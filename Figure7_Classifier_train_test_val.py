#!/usr/bin/env python3

import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings(action='ignore')

## Model-related imports
# Scaler
from sklearn.preprocessing import StandardScaler

# Model Selection
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import cross_validate

# Models
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier

# Metrics
from sklearn.metrics import accuracy_score
from sklearn.metrics import auc
from sklearn.metrics import average_precision_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import roc_curve

# Figure configurations
plt.style.use('ggplot')
plt.rcParams['grid.color'] = 'grey'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['xtick.color'] = 'black'
plt.rcParams['ytick.color'] = 'black'
plt.rcParams['axes.labelcolor'] = 'black'
# make fonts compatible with illustrator (prevents issues with vector graphics)
plt.rcParams['svg.fonttype'] = 'none'    


def import_train_data(train_data_file):
    """
    Import Training data set
    :param train_data_file <str>: Path to Training data set file in tsv format
    :return <pd.DataFrame>: Training data set
    """

    df = pd.read_table(train_data_file, index_col=0)
    columns = ['outcome', 'IL-1b', 'IL-6', 'TNF-a', 'CXCL9', 'CXCL10', 'CXCL11',
               'IFN-lambda', 'IL-8', 'IL-12p70', 'IFN-a', 'IL-28', 'GM-CSF', 
                'IFN-b', 'IL-10', 'IFN-g']

    df.columns = columns

    df.dropna(axis=0, inplace=True)

    # exclude Healthy Patient Samples from Dataset
    df = df[df['outcome'] != 'Healthy']

    # outcome values to bool
    df['outcome'] = df['outcome'] == 'Streptococcus pyogenes'
    
    return df


def import_validation_data(val_data_file):
    """
    Import Validation Data set
    :param val_data_file <str>: Path to Validation data set file in tsv format
    :return <pd.DataFrame>: Validation data set
    """

    df = pd.read_table(val_data_file, index_col=0)

    # data set contains only top 3 features
    columns = ['outcome', 'CXCL9', 'CXCL10', 'CXCL11']
    df.columns = columns

    df.dropna(axis=0, inplace=True)
    df['outcome'] = df['outcome'] == 'Streptococcus pyogenes'

    return df


def get_Xy(df_raw, outcome_feature, selected_features=None, standard_scale=False):
    """
    Get feature matrix and target vector from Dataframe
    :param df_raw <pd.DataFrame>:
    :param outcome_feature <str>: Name of target variable columns
    :param selected_features <list>: Specify to select a subset of features
    :param standard_scale <bool>: Standard scale data
    :return <tuple>: feature matrix <np.array>, target vector <np.array>
    """
    
    df = df_raw.copy()
    df.dropna(axis=0, inplace=True)
    
    y =  np.array(df.loc[:, outcome_feature])
    
    if selected_features:
        X = df.loc[:, selected_features]
    else:
        X = df.drop(outcome_feature, axis=1)
    
    if standard_scale == True:
        X = StandardScaler().fit_transform(X)
    
    X = np.array(X)
    
    return X, y


def get_cv_AUC_data_scores(X, y, classifier, cv=None, scores=None):
    """
    Evaluate Classifer performance by cross-validation
    :param X <np.array>: predictors
    :param y <np.1Darray>: dependent variable
    :param classifier <object>: initialized model that should be evaluated
    :param cv <cross validator object>: use specified cv or 5x2 Repeated stratified KFold cv
    :param scores <list>: list of sklearn scoring metrics
    :return: tuple containing cv-outcomes, cv-probabilities, pd.DataFrame of cv metrics
    """
    
    if not scores:
        scores = ['accuracy', 'precision', 'recall', 'roc_auc', 'f1','average_precision']
    
    if not cv:
        # 5x2 stratified CV
        cv = RepeatedStratifiedKFold(n_repeats=5, n_splits=2, random_state=random_state)
    
    cv_results = dict()
    y_real = []
    y_proba = []
    print('\tRunning Cross Validation')
    for i, (train, test) in enumerate(cv.split(X, y)):
        print('\t\tRound {}'.format(i+1))
        classifier.fit(X[train], y[train])
        pred_proba = classifier.predict_proba(X[test])
        y_real.append(y[test])
        y_proba.append(pred_proba[:,1])
        
        cv_scores = cross_validate(classifier, X=X, y=y, cv=cv, scoring=scores, return_train_score=False)
        for score_name in scores:
            # score names contain 'test_' prefix
            k = 'test_{}'.format(score_name)
            cv_results[score_name] = cv_scores[k]
        
    return (y_real, y_proba, pd.DataFrame(cv_results))
    

def export_cv_report(cv_scores, fname):
    """
    Export source data of cross validation experiment to Excel File
    :param cv_scores <dict>: Cross-Validation Scores
    :param fname <str>: Path to Excel File
    """
    
    writer = pd.ExcelWriter(fname, engine='xlsxwriter')
    stats_aggr = dict()

    for model_name, data in cv_scores.items():
        # export raw values
        data.to_excel(writer, sheet_name=model_name)
        
        # calculate statistics based on cv
        row_data = dict()
        row_data.update(data.mean(axis=0).rename(lambda x: '{}_mean'.format(x)))
        row_data.update(data.std(axis=0).rename(lambda x: '{}_std'.format(x)))
        row_data.update(data.sem(axis=0).rename(lambda x: '{}_sem'.format(x)))
        stats_aggr[model_name] = row_data
    
    df = pd.DataFrame.from_dict(stats_aggr, orient='index')
    df = df[sorted(df.columns)]
    df.to_excel(writer, sheet_name='statistics')
    writer.save()


def figure7A_evaluate_model(df_train_test, classifiers, cv, features, export_cv_data=True):
    """
    Compare RF, LR, and SVC Model Performance through Stratified Crossvalidation
    :param df_train_test <pd:DataFrame>: Train and Test data set
    :param classifiers <dict>: Classifiers to test {Name: initialized model}
    :param cv <cross validator object>: initialized cross validator
    :param export_cv_data <bool>: Export CV Validation metrics to Excel Table
    :return <tuple>: Scores of Cross-Validation runs for each tested model <dict>, CV run actuals and pred probabilities <dict>
    """

    # Analysis will be done on ALL features
    cv_scores = dict()
    AUC_raw_data = dict()

    for name, classifier in classifiers:
        print('\tRunning CV on {} classifier on all features'.format(name))

        # get cv metrics for AUROC
        if name.split(' ')[0] in ['SVC', 'Logistic Regression']:
            X, y = get_Xy(df_train_test, 'outcome', selected_features=features, standard_scale=True)
        else:
            X, y = get_Xy(df_train_test, 'outcome', selected_features=features, standard_scale=False)
            
        y_actuals, y_probas, df_scores = get_cv_AUC_data_scores(X, y, classifier, cv=cv)
        AUC_raw_data[name] = (y_actuals, y_probas)

        cv_scores[name] = df_scores

        print()
    
    if export_cv_data == True:
        print('Figure 7A source data exported')
        export_cv_report(cv_scores, './results/figure7/Fig7A_Model_comparison_ROC.xlsx')
    
    return cv_scores, AUC_raw_data


def figure7A_plot(cv_scores, AUC_raw_data):
    """
    Plot ROC from Cross Validation metrics
    :param cv_scores <dict>: Cross Validation metrics
    :param AUC_raw_data <dict>: actuals and prediction probabilities for each CV run and model
    """

    #fig, ax = plt.subplots(1, 1, figsize=(10,10))
    fig, ax = plt.subplots(1, 1)

    colors = {'Logistic Regression': '#e34b34',
              'SVC': '#1f77b5',
              'Random Forest': '#8dba41'}
    
    for i, name in enumerate(AUC_raw_data):
        y_real, y_proba = AUC_raw_data[name]
        y_real = np.array(y_real)
        y_proba = np.array(y_proba)
        
        ###########
        # ROC-AUC #
        ###########
        mean_fpr = np.linspace(0, 1, 100)
        tprs = []
        aucs = []

        # repeats x [fprs, tprs, thresholds]
        roc_curves = np.array([roc_curve(y_real[i], y_proba[i]) for i in range(y_real.shape[0])])

        for fpr, tpr, threshold in roc_curves:
            tprs.append(np.interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            roc_auc = auc(fpr, tpr)
            aucs.append(roc_auc)

        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        std_tpr = np.std(tprs, axis=0)

        mean_auc = auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)
        
        lab = '{} {:.3f} (+/- {:.2f})'.format(name, mean_auc, std_auc)
        ax.plot(mean_fpr, mean_tpr, label=lab, lw=2, c=colors[name])
    
    ###########
    # ROC-AUC #
    ###########
    ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
                label='Baseline', alpha=.8)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.legend(fontsize='small', bbox_to_anchor=(1.05, 1.05))
    ax.set_xlim([-0.05, 1.])
    ax.set_ylim([0., 1.05])
    ax.set_aspect('equal')
    ax.set_title('ROC curve')

    fig.tight_layout()
    print('Figure 7A created')
    fig.savefig('./results/figure7/Fig7A_Model_comparison_ROC.png', dpi=200)
    fig.savefig('./results/figure7/Fig7A_Model_comparison_ROC.svg', dpi=200)
    

def figure7C_evaluate_model(df_train_test, classifier, cv, feature_sets, export_cv_data=True):
    """
    Compare RF Model Performance on Full Panel or Sparse Data Set through Stratified Crossvalidation
    :param df_train_test <pd:DataFrame>: Train and Test data set
    :param classifier <object>: Initialized Classifier
    :param cv <cross validator object>: Initialized cross validator
    :param feature_sets <dict>: Features sets to test {Name: [features]}
    :param export_cv_data <bool>: Export CV Validation metrics to Excel Table
    :return <tuple>: Scores of Cross-Validation runs for each tested model <dict>, CV run actuals and pred probabilities <dict>
    """

    cv_scores = dict()
    AUC_raw_data = dict()

    for name, feature in feature_sets.items():
        print('\tRunning CV for Random Forest with {}'.format(name))
        # get data for AUROC and PR
        X, y = get_Xy(df_train_test, 'outcome', selected_features=feature, standard_scale=False)

        y_actuals, y_probas, df_scores = get_cv_AUC_data_scores(X, y, classifier, cv=cv)
        AUC_raw_data[name] = (y_actuals, y_probas)
        # get all metrics
        cv_scores[name] = df_scores
        print()
    
    if export_cv_data == True:
        export_cv_report(cv_scores, './results/figure7/Fig7C_Feat_Comparison_ROC_RF.xlsx')

    return cv_scores, AUC_raw_data


def figure7C_plot(cv_scores, AUC_raw_data):
    """
    Plot ROC from Cross Validation metrics
    :param cv_scores <dict>: Cross Validation metrics
    :param AUC_raw_data <dict>: actuals and prediction probabilities for each CV run and model
    """

    #fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    fig, ax = plt.subplots(1, 1)

    colors = {'Full Panel': '#e34b34',
              'Sparse Panel': '#8dba41'}

    for i, name in enumerate(AUC_raw_data):
        y_real, y_proba = AUC_raw_data[name]
        y_real = np.array(y_real)
        y_proba = np.array(y_proba)
        
        ###########
        # ROC-AUC #
        ###########
        mean_fpr = np.linspace(0, 1, 100)
        tprs = []
        aucs = []

        # repeats x [fprs, tprs, thresholds]
        roc_curves = np.array([roc_curve(y_real[i], y_proba[i]) for i in range(y_real.shape[0])])

        for fpr, tpr, threshold in roc_curves:
            tprs.append(np.interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            roc_auc = auc(fpr, tpr)
            aucs.append(roc_auc)

        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        std_tpr = np.std(tprs, axis=0)

        mean_auc = auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)
        
        lab = '{} {:.3f} (+/- {:.2f})'.format(name, mean_auc, std_auc)
        ax.plot(mean_fpr, mean_tpr, label=lab, lw=2, c=colors[name])


    ###########
    # ROC-AUC #
    ###########
    ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
                label='Baseline', alpha=.8)

    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.legend(fontsize='small', bbox_to_anchor=(1.05, 1.05))
    ax.set_xlim([-0.05, 1.])
    ax.set_ylim([0., 1.05])
    ax.set_aspect('equal')
    ax.set_title('ROC curve')

    fig.tight_layout()
    fig.savefig('./results/figure7/Fig7C_Feat_Comparison_ROC_RF.png', dpi=200)
    fig.savefig('./results/figure7/Fig7C_Feat_Comparison_ROC_RF.svg', dpi=200)


def figure7D_final_model(df_train_test, df_val, classifier, features):
    """
    :param df_train_test <pd:DataFrame>: Train and Test data set
    :param df_val <pd:DataFrame>: Validation data set (external cohort)
    :param classifier <object>: Initialized Classifier
    :param cv <cross validator object>: Initialized cross validator
    :param features <list>: Feature set to test
    """

    X_train, y_train = get_Xy(df_train_test, 'outcome', selected_features=features, standard_scale=False)
    X_val, y_val = get_Xy(df_val, 'outcome', selected_features=features, standard_scale=False)

    classifier.fit(X_train, y_train)

    y_proba = classifier.predict_proba(X_val)[:,1]
    y_pred = classifier.predict(X_val)

    # calculate metrics
    fpr, tpr, _ = roc_curve(y_val, y_proba)
    roc_auc = auc(fpr, tpr)
 
    #fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    fig, ax = plt.subplots(1, 1)
    lab = 'Final Model {:.3f}'.format(roc_auc)
    ax.plot(fpr, tpr, lw=2, c='#8dba41', label=lab)
    ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
                label='Baseline', alpha=.8)

    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.legend(fontsize='small', bbox_to_anchor=(1.05, 1.05))
    ax.set_xlim([-0.05, 1.])
    ax.set_ylim([0., 1.05])
    ax.set_aspect('equal')
    ax.set_title('ROC curve')

    fig.tight_layout()

    fig.savefig('./results/figure7/Fig7D_FinalModel_ROC.png', dpi=200)
    fig.savefig('./results/figure7/Fig7D_FinalModel_ROC.svg', dpi=200)

    print()
    print('Confusion Matrix:')
    print(confusion_matrix(y_val, y_pred))
    
    print('ROC-AUC:', auc(fpr, tpr))
    print('Average Precision:', average_precision_score(y_val, y_proba))
    print('Precision:', precision_score(y_val, y_pred))
    print('Recall:', recall_score(y_val, y_pred))
    print('F1:', f1_score(y_val, y_pred))
    print('Accuracy:', accuracy_score(y_val, y_pred))


def main():
    # global variables
    train_data_file = './data/fig7/train_data_all_markers.tsv'
    val_data_file = './data/fig7/external_cohort_data_top3-chemokines.tsv'

    random_state = 42
    repeats = 10

    # create output directory
    os.makedirs('results/figure7', exist_ok=True)

    best_params = {'LogisticRegression': {'C': 1, 'penalty': 'l2', 'solver': 'saga', 'tol': 1e-6},
                   'SVC': {'kernel': 'linear', 'C': 10, 'gamma': 1},
                   'RandomForest': {'max_features': 'sqrt', 'min_samples_leaf': 3, 'n_estimators': 250}}
    
    feature_sets = {'Full Panel': ['IL-1b', 'IL-6', 'TNF-a', 'CXCL9', 'CXCL10', 'CXCL11', 'IFN-lambda', 'IL-8', 'IL-12p70', 'IFN-a', 'IL-28', 'GM-CSF', 'IFN-b', 'IL-10', 'IFN-g'],
                    'Sparse Panel': ['CXCL10', 'CXCL9', 'CXCL11']}

    lr = LogisticRegression(**best_params['LogisticRegression'])
    svc = SVC(**best_params['SVC'], probability=True)
    rf = RandomForestClassifier(**best_params['RandomForest'], n_jobs=-1)
    
    classifiers = [('Logistic Regression', lr), ('SVC', svc), ('Random Forest', rf)]

    cv = RepeatedStratifiedKFold(n_repeats=repeats, n_splits=2, random_state=random_state)

    df_train = import_train_data(train_data_file)
    df_val = import_validation_data(val_data_file)

    #### FIGURE 7A ####
    print('#### Model Evaluation on Full data panel (Fig7A) ####')
    cv_scores, AUC_raw_data = figure7A_evaluate_model(df_train, classifiers, cv, feature_sets['Full Panel'], export_cv_data=True)
    figure7A_plot(cv_scores, AUC_raw_data)

    ####  FIGURE 7C ####
    print('#### Random Forest Evaluation on Full and Sparse data panel (Fig7C) ####')
    cv_scores, AUC_raw_data = figure7C_evaluate_model(df_train, rf, cv, feature_sets, export_cv_data=True)
    figure7C_plot(cv_scores, AUC_raw_data)

    #### Figure 7D ####
    print('#### Evaluation of Final Model on Full data panel (Fig7D) ####')
    rf = RandomForestClassifier(**best_params['RandomForest'], n_jobs=-1, random_state=random_state)
    figure7D_final_model(df_train, df_val, rf, feature_sets['Sparse Panel'])


if __name__ == '__main__':
    main()