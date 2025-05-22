
import os
import warnings
import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
from multiprocessing import cpu_count
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import (
    average_precision_score,
    roc_auc_score,
    accuracy_score,
    confusion_matrix,
    matthews_corrcoef
)
import parmap
import rpy2.robjects as robjects
from rpy2.robjects import r

# Suppress warnings
warnings.filterwarnings('ignore')

# Load R function
r.source('./code/PrOCTOR.R')
getStructuralFeatures = robjects.globalenv['getStructuralFeatures']

# Constants
CHEM_FEATURES = [
    'MolecularWeight', 'XLogP', 'HydrogenBondDonorCount', 'HydrogenBondAcceptorCount',
    'PolarSurfaceArea', 'FormalCharge', 'NumRings', 'RotatableBondCount', 'Refractivity',
    'Ro5', 'Ghose', 'Veber', 'wQED'
]
TARGET_FEATURES = [f"{x}_{y}" for x in ['Cell','Mouse','Human'] for y in ['Gene essentiality','Expression profile','Network connectivity']]

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', required=True, help='Evaluate model performance (e) or Predict new drug (p)')
    parser.add_argument('--input', help='Input data path for new drug prediction')
    parser.add_argument('--hg', default='hg38', help='Ensembl hg version of target genes')
    parser.add_argument('--output', required=True, help='Output directory path')
    return parser.parse_args()

def preprocess_data():
    data = pd.read_excel('./Supplementary_tables/Supplementary_Table_3.xlsx', header=1)
    data['Drug status'] = data['Drug status'].replace({'Risky drug': 1, 'Approved drug': 0}).astype(int)
    data = data.set_index('STITCH ID').sample(frac=1)
    data = data.rename(columns={
        'Gene essentiality': 'Cell_Gene essentiality',
        'Expression profile': 'Cell_Expression profile',
        'Network connectivity': 'Cell_Network connectivity',
        'Gene essentiality.1': 'Mouse_Gene essentiality',
        'Expression profile.1': 'Mouse_Expression profile',
        'Network connectivity.1': 'Mouse_Network connectivity',
        'Gene essentiality.2': 'Human_Gene essentiality',
        'Expression profile.2': 'Human_Expression profile',
        'Network connectivity.2': 'Human_Network connectivity'
    })
    return data, list(data.columns[:-1])

def evaluate_model(index, data, feature_cols):
    train_idx, test_idx = index
    X_train, X_test = data.loc[train_idx, feature_cols], data.loc[test_idx, feature_cols]
    y_train, y_test = data.loc[train_idx, 'Drug status'], data.loc[test_idx, 'Drug status']

    model = RandomForestClassifier(n_estimators=1000)
    model.fit(X_train, y_train)
    return model.predict_proba(X_test)[0, 1]

def evaluate_model_wrapper(index):
    args = parse_arguments()
    data, feature_cols = preprocess_data()
    return evaluate_model(index, data, feature_cols)

def evaluate_performance(data, feature_cols):
    args = parse_arguments()
    loo = LeaveOneOut()
    index_list = [[data.iloc[train].index.tolist(), data.iloc[test].index.tolist()]
                  for train, test in loo.split(data)]
    
    y_proba_list = parmap.map(evaluate_model_wrapper, index_list, pm_pbar=True, pm_processes=cpu_count())
    y_pred_list = [1 if x >= 0.4 else 0 for x in y_proba_list]
    y_true = data['Drug status'].tolist()

    metrics = {
        'AUPRC': average_precision_score(y_true, y_proba_list),
        'AUROC': roc_auc_score(y_true, y_proba_list),
        'Accuracy': accuracy_score(y_true, y_pred_list),
        'MCC': matthews_corrcoef(y_true, y_pred_list)
    }
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred_list).ravel()
    metrics.update({
        'Precision': tp / (tp + fp),
        'Recall': tp / (tp + fn),
        'Specificity': tn / (tn + fp)
    })

    os.makedirs(args.output, exist_ok=True)
    pd.DataFrame([metrics]).to_excel(os.path.join(args.output, 'evaluation_performance.xlsx'), index=False)
    return pd.DataFrame([metrics])

def load_target_matrices():
    mat_ch = pd.read_excel('./Supplementary_tables/Supplementary_Table_2.xlsx', sheet_name='Cell_Human', header=1)
    mat_ch = mat_ch.rename(columns={
        'Gene essentiality': 'Cell_Gene essentiality',
        'Expression profile': 'Cell_Expression profile',
        'Network connectivity': 'Cell_Network connectivity',
        'Gene essentiality.1': 'Human_Gene essentiality',
        'Expression profile.1': 'Human_Expression profile',
        'Network connectivity.1': 'Human_Network connectivity'
    }).dropna().reset_index(drop=True)

    mat_m = pd.read_excel('./Supplementary_tables/Supplementary_Table_2.xlsx', sheet_name='Mouse', header=1)
    mat_m = mat_m.rename(columns={
        'Gene essentiality': 'Mouse_Gene essentiality',
        'Expression profile': 'Mouse_Expression profile',
        'Network connectivity': 'Mouse_Network connectivity'
    }).dropna().reset_index(drop=True)

    return mat_ch, mat_m

def extract_target_features(entry, mat_ch, mat_m, hg_version):
    h_targets = entry[2].split('|')
    m_targets = entry[3].split('|')
    mat_ch_temp = mat_ch[mat_ch[f'{hg_version}_gene_id'].isin(h_targets)]
    mat_m_temp = mat_m[mat_m['mouse_gene_id'].isin(m_targets)]
    
    feature_dict = {}
    for feat in TARGET_FEATURES:
        if 'Mouse' in feat:
            feature_dict[feat] = mat_m_temp[feat].max()
        else:
            feature_dict[feat] = mat_ch_temp[feat].max()

    feature_dict['Human_target_coverage'] = 'PASS' if len(mat_ch_temp) / len(h_targets) >= 0.5 else 'FAIL'
    feature_dict['Mouse_target_coverage'] = 'PASS' if len(mat_m_temp) / len(m_targets) >= 0.5 else 'FAIL'

    return feature_dict

def extract_target_features_wrapper(row):
    args = parse_arguments()
    mat_ch, mat_m = load_target_matrices()
    return extract_target_features(row, mat_ch, mat_m, args.hg)

def extract_structural_features(smiles):
    raw = getStructuralFeatures(SMILE=smiles)
    feats = list(raw)
    feature_dict = {'SMILES_featurization': 'PASS'}

    for idx, feat in enumerate(CHEM_FEATURES):
        val = feats[idx]
        feature_dict[feat] = val
        if np.isnan(val):
            feature_dict['SMILES_featurization'] = 'FAIL'
    
    return feature_dict

def predict_new_drugs(args, data, feature_cols):
    new_data = pd.read_csv(args.input, sep='\t')

    print('Featurizing SMILES...')
    structural_features = [extract_structural_features(row[1]) for row in tqdm(new_data.to_numpy())]
    target_features = parmap.map(extract_target_features_wrapper, new_data.to_numpy(), pm_pbar=True, pm_processes=cpu_count())

    for i in range(len(structural_features)):
        target_features[i].update(structural_features[i])
        target_features[i]['Drug'] = new_data['Drug'].iloc[i]

    df_features = pd.DataFrame(target_features)
    filter_pass = (
        (df_features['Human_target_coverage'] == 'PASS') &
        (df_features['Mouse_target_coverage'] == 'PASS') &
        (df_features['SMILES_featurization'] == 'PASS')
    )
    df_filtered = df_features[filter_pass].copy()
    df_failed = df_features[~filter_pass].copy()

    model = RandomForestClassifier(n_estimators=1000)
    model.fit(data[feature_cols], data['Drug status'])

    X_new = df_filtered[feature_cols]
    y_proba = model.predict_proba(X_new)[:, 1]
    y_pred = [1 if p >= 0.4 else 0 for p in y_proba]

    df_filtered['Toxicity probability'] = y_proba
    df_filtered['Predicted drug status'] = ['Risky drug' if p == 1 else 'Approved drug' for p in y_pred]

    result = df_filtered[['Drug', 'Toxicity probability', 'Predicted drug status']]
    for col in ['Human_target_coverage', 'Mouse_target_coverage', 'SMILES_featurization']:
        result[col.replace('_', ' ')] = 'PASS'

    if not df_failed.empty:
        df_failed_out = df_failed[['Drug', 'Human_target_coverage', 'Mouse_target_coverage', 'SMILES_featurization']]
        df_failed_out.columns = ['Drug', 'Human target coverage', 'Mouse target coverage', 'SMILES featurization']
        result = pd.concat([result, df_failed_out], ignore_index=True)

    os.makedirs(args.output, exist_ok=True)
    result.to_excel(os.path.join(args.output, f"{os.path.basename(args.input)}.xlsx"), index=False)

def main():
    args = parse_arguments()
    data, feature_cols = preprocess_data()

    if args.mode == 'e':
        perf_df = evaluate_performance(data, feature_cols)
        print(perf_df)
    elif args.mode == 'p':
        predict_new_drugs(args, data, feature_cols)

if __name__ == '__main__':
    main()
