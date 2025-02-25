# Import required libraries
import os
import pandas as pd
import numpy as np
import xgboost as xgb
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score, accuracy_score, roc_curve, precision_recall_curve, auc
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer

# Set working directory and read data
os.chdir("Your own path")
data = pd.read_csv("Your own CSV")
target_column = 'Group'
data[target_column] = data[target_column].astype('category')

# Split features and target
X = data.drop(columns=target_column)
y = data[target_column]

# Define preprocessor for numerical features
preprocessor = ColumnTransformer(
    transformers=[
        ('num', SimpleImputer(strategy='median'), X.columns)
    ]
)

# Define dictionary of models to be evaluated
models = {
    'Logistic Regression': LogisticRegression(max_iter=10000, solver='liblinear'),
    'Random Forest': RandomForestClassifier(random_state=121),
    'AdaBoost': AdaBoostClassifier(random_state=121),
    'GradientBoosting': GradientBoostingClassifier(random_state=121),
    'XGBoost': xgb.XGBClassifier(use_label_encoder=False, eval_metric='logloss', seed=121),
    'KNN': KNeighborsClassifier(),
    'MLP': MLPClassifier(max_iter=2000, random_state=121),
    'Decision Tree': DecisionTreeClassifier(random_state=121),
    'SVM': SVC(probability=True, random_state=121)
}

# Set color cycle for plotting
colors = plt.cm.Set3(np.linspace(0, 1, len(models)))

# Initialize dictionaries to store results and predictions
results = {}
predictions = {}

# Train models and collect predictions
for (name, model), color in zip(models.items(), colors):
    # Create pipeline with preprocessing and model
    pipeline = Pipeline(steps=[('preprocessor', preprocessor),
                             ('scaler', StandardScaler()),
                             ('classifier', model)])
    
    # Fit model and get predictions
    pipeline.fit(X, y)
    y_pred = pipeline.predict_proba(X)[:, 1]
    predictions[name] = y_pred
    
    # Calculate performance metrics
    results[name] = {
        'auc': roc_auc_score(y, y_pred),
        'acc': accuracy_score(y, pipeline.predict(X))
    }

# Plot ROC curve
plt.figure(figsize=(10, 8))
for (name, y_pred), color in zip(predictions.items(), colors):
    fpr, tpr, _ = roc_curve(y, y_pred)
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, color=color, label=f'{name} (AUC = {roc_auc:.2f})')

plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve')
plt.legend(loc='lower right')
plt.grid(True)
plt.show()

# Plot Precision-Recall curve
plt.figure(figsize=(10, 8))
for (name, y_pred), color in zip(predictions.items(), colors):
    precision, recall, _ = precision_recall_curve(y, y_pred)
    prc_auc = auc(recall, precision)
    plt.plot(recall, precision, color=color, label=f'{name} (AUC = {prc_auc:.2f})')

plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve')
plt.legend(loc='lower left')
plt.grid(True)
plt.show()

# Define function for Decision Curve Analysis
def calculate_net_benefit(y_true, y_pred_proba, threshold):
    """
    Calculate net benefit for a given threshold
    Args:
        y_true: True labels
        y_pred_proba: Predicted probabilities
        threshold: Decision threshold
    Returns:
        Net benefit value
    """
    y_pred = (y_pred_proba >= threshold).astype(int)
    TP = np.sum((y_true == 1) & (y_pred == 1))
    FP = np.sum((y_true == 0) & (y_pred == 1))
    n = len(y_true)
    
    if TP + FP == 0:
        return 0
    
    net_benefit = (TP/n) - (FP/n) * (threshold/(1-threshold))
    return net_benefit

def plot_dca_curve(y_true, predictions_dict, thresholds=np.arange(0, 1.01, 0.01),
                  title="Decision Curve Analysis", colors=None):
    """
    Plot Decision Curve Analysis
    Args:
        y_true: True labels
        predictions_dict: Dictionary of model predictions
        thresholds: Array of threshold values
        title: Plot title
        colors: Color scheme for plotting
    """
    y_true_num = pd.get_dummies(y_true).iloc[:, 1].values
    
    plt.figure(figsize=(10, 8))
    
    if colors is None:
        colors = plt.cm.Set3(np.linspace(0, 1, len(predictions_dict)))
    
    # Calculate and plot "treat all" net benefit
    all_treat = [np.mean(y_true_num) - threshold/(1-threshold)*(1-np.mean(y_true_num))
                 for threshold in thresholds]
    plt.plot(thresholds, all_treat, 'k--', label='Treat All')
    
    # Plot "treat none" baseline
    plt.plot(thresholds, np.zeros_like(thresholds), 'k-', label='Treat None')
    
    # Calculate and plot net benefit for each model
    for (name, y_pred), color in zip(predictions_dict.items(), colors):
        net_benefits = []
        for threshold in thresholds:
            nb = calculate_net_benefit(y_true_num, y_pred, threshold)
            net_benefits.append(nb)
        plt.plot(thresholds, net_benefits, color=color, label=name)
    
    plt.xlim(0, 1)
    plt.ylim(-0.05, max(all_treat) + 0.1)
    plt.xlabel('Threshold Probability')
    plt.ylabel('Net Benefit')
    plt.title(title)
    plt.legend(loc='lower left')
    plt.grid(True)
    
    return plt.gcf()

# Plot DCA curve
dca = plot_dca_curve(y, predictions)
plt.show()

# Calculate and display comprehensive metrics
from sklearn.metrics import precision_score, recall_score, f1_score

metrics = {}
for name, model in models.items():
    pipeline = Pipeline(steps=[('preprocessor', preprocessor),
                             ('scaler', StandardScaler()),
                             ('classifier', model)])
    
    pipeline.fit(X, y)
    y_pred_class = pipeline.predict(X)
    y_pred_prob = pipeline.predict_proba(X)[:, 1]
    
    metrics[name] = {
        'Accuracy': accuracy_score(y, y_pred_class),
        'Precision': precision_score(y, y_pred_class),
        'Recall': recall_score(y, y_pred_class),
        'F1-score': f1_score(y, y_pred_class),
        'ROC-AUC': roc_auc_score(y, y_pred_prob)
    }

# Create and display metrics DataFrame
metrics_df = pd.DataFrame(metrics).T
print("\nModel Performance Metrics:")
print(metrics_df.round(3))

# Save results to CSV
metrics_df.to_csv('model_metrics.csv')

# Plot radar charts for each metric
metrics_list = ['Accuracy', 'Precision', 'Recall', 'F1-score', 'ROC-AUC']
model_names = list(models.keys())
angles = np.linspace(0, 2*np.pi, len(model_names), endpoint=False)
angles = np.concatenate((angles, [angles[0]]))
model_names_plot = np.concatenate((model_names, [model_names[0]]))

def add_value_labels(ax, angles, values):
    """
    Add value labels to radar chart
    """
    for angle, value in zip(angles[:-1], values[:-1]):
        ha = 'left' if 0 <= angle <= np.pi else 'right'
        offset = 0.1 if 0 <= angle <= np.pi else -0.1
        ax.text(angle, value + 0.05, f'{value:.3f}',
                ha=ha, va='center')

# Plot radar chart for each metric
for metric in metrics_list:
    plt.figure(figsize=(10, 8))
    ax = plt.subplot(111, projection='polar')
    
    values = [metrics[model][metric] for model in model_names]
    values = np.concatenate((values, [values[0]]))
    
    ax.plot(angles, values, 'o-', linewidth=2, label=metric)
    ax.fill(angles, values, alpha=0.25)
    
    add_value_labels(ax, angles, values)
    
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(model_names, fontsize=10)
    ax.set_ylim(0, 1)
    ax.set_title(f'{metric} Performance')
    ax.grid(True)
    
    plt.tight_layout()
    plt.show()
