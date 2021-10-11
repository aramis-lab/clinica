import numpy as np
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    classification_report,
)


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

    accuracy = (true_positive + true_negative) / (
        true_positive + true_negative + false_positive + false_negative
    )

    if (true_positive + false_negative) != 0:
        sensitivity = true_positive / (true_positive + false_negative)
    else:
        sensitivity = 0.0

    if (false_positive + true_negative) != 0:
        specificity = true_negative / (false_positive + true_negative)
    else:
        specificity = 0.0

    if (true_positive + false_positive) != 0:
        ppv = true_positive / (true_positive + false_positive)
    else:
        ppv = 0.0

    if (true_negative + false_negative) != 0:
        npv = true_negative / (true_negative + false_negative)
    else:
        npv = 0.0

    balanced_accuracy = (sensitivity + specificity) / 2

    results = {
        "accuracy": accuracy,
        "balanced_accuracy": balanced_accuracy,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "ppv": ppv,
        "npv": npv,
        "confusion_matrix": {
            "tp": len(tp),
            "tn": len(tn),
            "fp": len(fp),
            "fn": len(fn),
        },
    }

    return results


def gram_matrix_linear(data):
    return np.dot(data, data.transpose())


def evaluate_prediction_multiclass(y, y_hat):

    balanced_accuracy = balanced_accuracy_score(y, y_hat)
    accuracy = accuracy_score(y, y_hat)

    results = {"accuracy": accuracy, "balanced_accuracy": balanced_accuracy}

    return results


def metric_distribution(
    metric, labels, output_path, num_classes=2, metric_label="balanced accuracy"
):
    """

    Distribution plots of various metrics such as balanced accuracy!

    metric is expected to be ndarray of size [num_repetitions, num_datasets]

    """
    # from __future__ import print_function, division

    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import cm
    from matplotlib.backends.backend_pdf import PdfPages

    num_repetitions = metric.shape[0]
    num_datasets = metric.shape[1]
    assert len(labels) == num_datasets, "Differing number of features and labels!"
    method_ticks = 1.0 + np.arange(num_datasets)

    fig, ax = plt.subplots(figsize=[9, 9])
    line_coll = ax.violinplot(
        metric,
        widths=0.8,
        bw_method=0.2,
        showmedians=True,
        showextrema=False,
        positions=method_ticks,
    )

    cmap = cm.get_cmap("Paired", num_datasets)
    for cc, ln in enumerate(line_coll["bodies"]):
        ln.set_facecolor(cmap(cc))
        ln.set_label(labels[cc])

    plt.legend(loc=2, ncol=num_datasets)

    ax.tick_params(axis="both", which="major", labelsize=15)
    ax.grid(axis="y", which="major")

    lower_lim = np.round(np.min([np.float64(0.9 / num_classes), metric.min()]), 3)
    upper_lim = np.round(np.max([1.01, metric.max()]), 3)
    step_tick = 0.1
    ax.set_ylim(lower_lim, upper_lim)

    ax.set_xticks(method_ticks)
    ax.set_xlim(np.min(method_ticks) - 1, np.max(method_ticks) + 1)
    ax.set_xticklabels(labels, rotation=45)  # 'vertical'

    ax.set_yticks(np.arange(lower_lim, upper_lim, step_tick))
    ax.set_yticklabels(np.arange(lower_lim, upper_lim, step_tick))
    # plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(metric_label, fontsize=16)

    fig.tight_layout()

    pp1 = PdfPages(output_path + ".pdf")
    pp1.savefig()
    pp1.close()

    return
