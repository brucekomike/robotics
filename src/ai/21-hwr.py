import tkinter as tk
from tkinter import filedialog
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import math

def load_digit_data(filepath):
    """
    加载数字数据。假设第一列是标签，其余是像素特征。
    """
    if not filepath:
        print("未选择文件。正在退出。")
        return None, None, None

    try:
        # 假设没有表头，数据是数字。
        # 使用 float32 以提高像素数据的内存效率。
        data = pd.read_csv(filepath, header=None, dtype=np.float32, delim_whitespace=True)
        print(f"数据加载成功。形状: {data.shape}")
    except Exception as e:
        print(f"从 {filepath} 加载数据时出错: {e}")
        return None, None, None

    if data.shape[1] < 2:
        print("错误：数据必须至少有两列（标签 + 一个像素特征）。")
        return None, None, None

    # 第一列是标签，其余是特征
    y = data.iloc[:, 0].astype(int).values
    X = data.iloc[:, 1:].values

    # 推断图像尺寸（假设是方形图像）
    num_features = X.shape[1]
    image_side = int(math.sqrt(num_features))
    if image_side * image_side != num_features:
        print(f"警告：特征数量 ({num_features}) 不是一个完全平方数。无法重塑为方形图像进行显示。")
        image_side = None # 无法确定非方形图像的边长
    else:
        print(f"推断的图像尺寸: {image_side}x{image_side}")

    return X, y, image_side

def plot_sample_digits(X, y, image_side, n_samples=5):
    """
    从数据集中绘制一些样本数字。
    """
    if image_side is None:
        print("由于无法确定图像尺寸或图像非方形，无法绘制样本数字。")
        return

    plt.figure(figsize=(10, 4))
    for i in range(n_samples):
        if i >= len(X):
            break
        plt.subplot(1, n_samples, i + 1)
        try:
            image = X[i].reshape(image_side, image_side)
            plt.imshow(image, cmap='gray')
            plt.title(f"标签: {y[i]}")
            plt.axis('off')
        except ValueError as e:
            print(f"重塑图像 {i} 以进行绘制时出错: {e}。特征数: {X[i].shape[0]}, 期望: {image_side*image_side}")
            plt.text(0.5, 0.5, '错误\n重塑中', ha='center', va='center')
            plt.axis('off')

    plt.suptitle(f"样本数字 ({image_side}x{image_side})")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # 调整布局为标题腾出空间
    plt.show()

def train_and_evaluate_classifier(X_train, y_train, X_test, y_test, hidden_layer_sizes, activation, model_name, solver='adam', max_iter=300, random_state=42):
    """
    训练一个 MLPClassifier 并打印评估指标。
    返回训练好的模型及其损失曲线。
    """
    print(f"\n--- 训练: {model_name} ---")
    print(f"激活函数: {activation}, 隐藏层: {hidden_layer_sizes}")

    mlp = MLPClassifier(
        hidden_layer_sizes=hidden_layer_sizes,
        activation=activation,
        solver=solver,
        max_iter=max_iter,
        random_state=random_state,
        learning_rate_init=0.001,
        early_stopping=True, # 启用早停
        n_iter_no_change=10, # 如果10次迭代后没有改善则停止
        tol=1e-4
    )

    mlp.fit(X_train, y_train)

    y_pred_train = mlp.predict(X_train)
    y_pred_test = mlp.predict(X_test)

    train_accuracy = accuracy_score(y_train, y_pred_train)
    test_accuracy = accuracy_score(y_test, y_pred_test)

    print(f"  训练集准确率: {train_accuracy:.4f}")
    print(f"  测试集准确率: {test_accuracy:.4f}")
    # print("\n  测试集分类报告:")
    # print(classification_report(y_test, y_pred_test, zero_division=0))
    # print("\n  测试集混淆矩阵:")
    # print(confusion_matrix(y_test, y_pred_test))

    return mlp, train_accuracy, test_accuracy, mlp.loss_curve_ if hasattr(mlp, 'loss_curve_') else None

def main_classification_task():
    """
    手写数字分类任务的主函数。
    """
    root = tk.Tk()
    root.withdraw() # 隐藏主 tkinter 窗口

    filepath = filedialog.askopenfilename(
        title="选择手写数字数据文件",
        filetypes=(("CSV files", "*.csv"), ("Text files", "*.txt"), ("Data files", "*.dat"), ("All files", "*.*"))
    )

    if not filepath:
        print("未选择文件。正在退出。")
        return

    X, y, image_side = load_digit_data(filepath)

    if X is None or y is None:
        return

    # 显示一些样本数字
    if image_side is not None:
        plot_sample_digits(X[:20], y[:20], image_side, n_samples=5) # 绘制前20个样本中的前5个

    # 预处理：缩放特征
    # 像素值通常在 0-255 或 0-1 之间。StandardScaler 进行中心化和缩放。
    # 如果你的数据已经是 0-1，缩放可能不是严格必要的，但通常不会有害。
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    print(f"\n数据已缩放。X_scaled 形状: {X_scaled.shape}")

    # 分割数据
    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.25, random_state=42, stratify=y)
    print(f"数据已分割: 训练集大小: {X_train.shape[0]}, 测试集大小: {X_test.shape[0]}")

    # --- 任务1：比较固定隐藏层大小下的 ReLU 和 Logistic ---
    # (例如，一个包含100个神经元的隐藏层，根据OCR第12页“包含100个隐藏节点的三层感知机网络”)
    # 这里的“三层”指的是输入层、一个隐藏层、输出层。
    print("\n\n=== 任务1：比较 ReLU 和 Logistic (100个隐藏节点) ===")
    fixed_hidden_nodes = (100,)
    
    model_relu_fixed, _, _, loss_relu_fixed = train_and_evaluate_classifier(
        X_train, y_train, X_test, y_test,
        hidden_layer_sizes=fixed_hidden_nodes,
        activation='relu',
        model_name="ReLU (100 个节点)"
    )

    model_logistic_fixed, _, _, loss_logistic_fixed = train_and_evaluate_classifier(
        X_train, y_train, X_test, y_test,
        hidden_layer_sizes=fixed_hidden_nodes,
        activation='logistic',
        model_name="Logistic (100 个节点)"
    )

    # 绘制固定隐藏层比较的损失曲线
    if loss_relu_fixed is not None and loss_logistic_fixed is not None:
        plt.figure(figsize=(10, 6))
        plt.plot(loss_relu_fixed, label='ReLU 损失 (100 个节点)')
        plt.plot(loss_logistic_fixed, label='Logistic 损失 (100 个节点)')
        plt.title('训练损失曲线 (100个隐藏节点)')
        plt.xlabel('迭代次数 (Epochs)')
        plt.ylabel('损失 (Loss)')
        plt.legend()
        plt.grid(True)
        plt.show()

    # --- 任务2：绘制训练/测试错误率与隐藏节点数的关系图 ---
    # (OCR 第12页：“隐藏节点数(由1-19)”)
    # 我们将使用单个隐藏层并改变其大小。
    # 错误率 = 1 - 准确率
    print("\n\n=== 任务2：错误率 vs. 隐藏节点数 (单层) ===")
    
    # 对于数字识别，1-19 的范围可能太小，无法获得良好性能。
    # 如果结果不佳，可以考虑一个稍大的范围，例如，以步长方式增加到50或100。
    # 目前，按照 OCR 的要求使用 1-19，但增加一些点以获得更好的曲线。
    hidden_node_counts = list(range(1, 20)) + [25, 30, 40, 50] # 扩展范围
    # hidden_node_counts = list(range(1, 20)) # 严格按照 OCR 的 1-19

    train_errors_relu = []
    test_errors_relu = []
    train_errors_logistic = []
    test_errors_logistic = []

    for n_nodes in hidden_node_counts:
        print(f"\n-- 使用 {n_nodes} 个隐藏节点进行测试 --")
        current_hidden_layer = (n_nodes,)

        # ReLU 激活函数
        _, train_acc_relu, test_acc_relu, _ = train_and_evaluate_classifier(
            X_train, y_train, X_test, y_test,
            hidden_layer_sizes=current_hidden_layer,
            activation='relu',
            model_name=f"ReLU ({n_nodes} 个节点)",
            max_iter=200 # 减少 max_iter 以加快循环速度，如果需要可以调整
        )
        train_errors_relu.append(1 - train_acc_relu)
        test_errors_relu.append(1 - test_acc_relu)

        # Logistic 激活函数
        _, train_acc_logistic, test_acc_logistic, _ = train_and_evaluate_classifier(
            X_train, y_train, X_test, y_test,
            hidden_layer_sizes=current_hidden_layer,
            activation='logistic',
            model_name=f"Logistic ({n_nodes} 个节点)",
            max_iter=200 # 减少 max_iter
        )
        train_errors_logistic.append(1 - train_acc_logistic)
        test_errors_logistic.append(1 - test_acc_logistic)

    # 绘制错误率曲线
    plt.figure(figsize=(14, 8))

    plt.subplot(1, 2, 1)
    plt.plot(hidden_node_counts, train_errors_relu, marker='o', linestyle='-', label='ReLU 训练错误率')
    plt.plot(hidden_node_counts, test_errors_relu, marker='x', linestyle='--', label='ReLU 测试错误率')
    plt.title('ReLU: 错误率 vs. 隐藏节点数')
    plt.xlabel('隐藏节点数量 (单层)')
    plt.ylabel('错误率 (1 - 准确率)')
    plt.legend()
    plt.grid(True)
    plt.xticks(hidden_node_counts[::2]) # 为了清晰起见，每隔一个计数显示一个刻度

    plt.subplot(1, 2, 2)
    plt.plot(hidden_node_counts, train_errors_logistic, marker='o', linestyle='-', label='Logistic 训练错误率')
    plt.plot(hidden_node_counts, test_errors_logistic, marker='x', linestyle='--', label='Logistic 测试错误率')
    plt.title('Logistic: 错误率 vs. 隐藏节点数')
    plt.xlabel('隐藏节点数量 (单层)')
    plt.ylabel('错误率 (1 - 准确率)')
    plt.legend()
    plt.grid(True)
    plt.xticks(hidden_node_counts[::2])

    plt.tight_layout()
    plt.show()

    print("\n\n手写数字识别任务完成。")

if __name__ == "__main__":
    main_classification_task()