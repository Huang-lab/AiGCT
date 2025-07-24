import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, f1_score

# Example: y_true and y_scores from your classifier
# y_true: true binary labels
# y_scores: model's predicted probability for the positive class
y_true = np.array([0, 1, 1, 0, 1, 0, 1, 0, 1, 1])
y_scores = np.array([0.1, 0.9, 0.8, 0.4, 0.95, 0.3, 0.6, 0.2, 0.85, 0.7])

# Get precision, recall, and thresholds
precision, recall, thresholds = precision_recall_curve(y_true, y_scores)

# Compute F1 scores
f1_scores = 2 * (precision * recall) / (precision + recall + 1e-8)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(thresholds, precision[:-1], label='Precision', color='b')
plt.plot(thresholds, recall[:-1], label='Recall', color='g')
plt.plot(thresholds, f1_scores[:-1], label='F1 Score', color='r')
plt.xlabel('Decision Threshold')
plt.ylabel('Score')
plt.title('Performance Metrics vs. Decision Threshold')
plt.legend()
plt.grid(True)
plt.show()
