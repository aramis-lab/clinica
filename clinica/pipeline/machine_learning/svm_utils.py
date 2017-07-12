
import numpy as np
import pandas as pd
import csv
from ntpath import basename, splitext
from scipy.spatial.distance import squareform
from sklearn.metrics import roc_auc_score


def evaluate_prediction(y, y_hat):

    true_positive = 0.0
    true_negative = 0.0
    false_positive = 0.0
    false_negative = 0.0

    tp = []
    tn = []
    fp = []
    fn = []

    for i in range(len(y)):
        if y[i] == 1:
            if y_hat[i] == 1:
                true_positive += 1
                tp.append(i)
            else:
                false_negative += 1
                fn.append(i)
        else:  # -1
            if y_hat[i] == 0:
                true_negative += 1
                tn.append(i)
            else:
                false_positive += 1
                fp.append(i)

    accuracy = (true_positive + true_negative) / (true_positive + true_negative + false_positive + false_negative)
    sensitivity = true_positive / (true_positive + false_negative)
    specificity = true_negative / (false_positive + true_negative)

    if (true_positive + false_positive) !=  0:
        ppv = true_positive / (true_positive + false_positive)
    else:
        ppv = 0

    if (true_negative + false_negative) != 0:
        npv = true_negative / (true_negative + false_negative)
    else:
        npv = 0

    balanced_accuracy = (sensitivity + specificity) / 2

    results = {'accuracy': accuracy,
               'balanced_accuracy': balanced_accuracy,
               'sensitivity': sensitivity,
               'specificity': specificity,
               'ppv': ppv,
               'npv': npv,
               'predictions': (tp, tn, fp, fn)
               }

    return results


def gram_matrix_linear(data):
    return np.dot(data, data.transpose())


def calculate_auc(svc, shared_x, train_indices, test_indices, y_test):

    x_train = shared_x[train_indices[np.array(svc.support_)], :]
    weights = np.sum(x_train * svc.dual_coef_.transpose(), 0)

    x_test = shared_x[test_indices, :]
    y_hat = x_test.dot(weights) + svc.intercept_

    return roc_auc_score(y_test, y_hat)


def save_subjects_prediction(subjects, diagnosis, y, y_hat, output_file):

    with open(output_file, 'w') as csvfile:

        fieldnames = ['Subject', 'Diagnose', 'Class_label', 'Predicted_label', 'Correct']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for i in range(len(subjects)):
            s = splitext(basename(subjects[i]))[0]
            dx = diagnosis[i]
            writer.writerow({'Subject': s,
                             'Diagnose': dx,
                             'Class_label': y[i],
                             'Predicted_label': int(y_hat[i]),
                             'Correct': int(y[i] == y_hat[i])
                             })


def results_to_tsv(results, diagnose_list, output_file):

    with open(output_file, 'w') as output:

        s = 'AUC\n'
        auc_list = [round(res[1]['auc'], 2) for res in sorted(results.items())]
        df = pd.DataFrame(squareform(auc_list), index=diagnose_list, columns=diagnose_list)
        print df
        s += df.to_csv(sep='\t')

        s += '\nBalanced Accuracy\n'
        balanced_accuracy_list = [round(res[1]['balanced_accuracy'], 2) for res in sorted(results.items())]
        df = pd.DataFrame(squareform(balanced_accuracy_list), index=diagnose_list, columns=diagnose_list)
        print df
        s += df.to_csv(sep='\t')

        s += '\nAccuracy\n'
        accuracy_list = [round(res[1]['accuracy'], 2) for res in sorted(results.items())]
        df = pd.DataFrame(squareform(accuracy_list), index=diagnose_list, columns=diagnose_list)
        print df
        s += df.to_csv(sep='\t')

        s += '\nSensitivity\n'
        sensitivity_list = [round(res[1]['sensitivity'], 2) for res in sorted(results.items())]
        df = pd.DataFrame(squareform(sensitivity_list), index=diagnose_list, columns=diagnose_list)
        print df
        s += df.to_csv(sep='\t')

        s += '\nSpecificity\n'
        specificity_list = [round(res[1]['specificity'], 2) for res in sorted(results.items())]
        df = pd.DataFrame(squareform(specificity_list), index=diagnose_list, columns=diagnose_list)
        print df
        s += df.to_csv(sep='\t')

        output.write(s)
