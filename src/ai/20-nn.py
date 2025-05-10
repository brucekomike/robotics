import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog

def load_and_preprocess_data(file_path):
    """
    加载CSV数据，进行预处理，选择特征和目标变量。
    """
    if not file_path: # 如果没有选择文件
        print("没有选择文件。程序将退出。")
        return None, None, None, None
        
    try:
        data = pd.read_csv(file_path, encoding='utf-8')
    except UnicodeDecodeError:
        try:
            data = pd.read_csv(file_path, encoding='gbk')
        except Exception as e:
            print(f"错误: 无法使用UTF-8或GBK编码读取文件 '{file_path}'。请检查文件编码。错误详情: {e}")
            return None, None, None, None
    except FileNotFoundError:
        print(f"错误: 文件 '{file_path}' 未找到。这不应该发生，因为文件是通过对话框选择的。")
        return None, None, None, None
    except Exception as e:
        print(f"读取文件 '{file_path}' 时发生未知错误: {e}")
        return None, None, None, None


    print("原始数据列名:", data.columns.tolist())

    # 确定PM2.5列名
    pm25_col_name = None
    possible_pm25_names = ['PM2.5', 'PM25', 'pm2.5', 'pm25'] # 常见PM2.5列名变体
    for name in possible_pm25_names:
        if name in data.columns:
            pm25_col_name = name
            break
    
    if pm25_col_name is None:
        print(f"错误：在CSV文件中未找到 PM2.5 列 (尝试过: {possible_pm25_names})。请检查列名。")
        return None, None, None, None

    target_col = pm25_col_name
    
    # 定义期望的特征列
    expected_feature_cols = ['AQI', 'PM10', 'SO2', 'CO', 'NO2', 'O3']
    actual_feature_cols = []
    for col in expected_feature_cols:
        if col in data.columns:
            if data[col].dtype in [np.number, np.int64, np.float64]: # 确保是数值类型
                 actual_feature_cols.append(col)
            else: # 尝试转换为数值类型
                try:
                    data[col] = pd.to_numeric(data[col], errors='raise')
                    actual_feature_cols.append(col)
                except (ValueError, TypeError):
                    print(f"警告：特征列 '{col}' 不是数值类型且无法转换，将被忽略。")
        else:
            print(f"警告：期望的特征列 '{col}' 未在数据中找到。")
    
    if not actual_feature_cols:
        print("错误：没有找到可用的数值特征列。")
        return None, None, None, None

    print(f"使用的特征列: {actual_feature_cols}")
    print(f"使用的目标列: {target_col}")

    X_df = data[actual_feature_cols].copy()
    y_series = data[target_col].copy()

    # 将目标列转换为数值类型，处理可能存在的非数值（例如占位符）
    y_series = pd.to_numeric(y_series, errors='coerce')

    # 处理特征和目标中的缺失值
    # 首先移除目标列为NaN的行，因为这些行无法用于训练或评估
    valid_indices = ~y_series.isna()
    y_series = y_series[valid_indices]
    X_df = X_df[valid_indices]

    if y_series.empty:
        print("错误：在目标列PM2.5中没有有效的数值数据。")
        return None, None, None, None

    # 对特征列填充缺失值
    for col in X_df.columns:
        if X_df[col].isnull().any():
            print(f"特征列 '{col}' 中存在缺失值，将使用均值填充。")
            X_df[col] = X_df[col].fillna(X_df[col].mean())

    # 数据标准化
    scaler_X = StandardScaler()
    X_scaled = scaler_X.fit_transform(X_df)
    
    return X_scaled, y_series.values, actual_feature_cols, scaler_X

def train_evaluate_plot_model(X_train, X_test, y_train, y_test, model_name, mlp_model, plot_loss=True, plot_predictions=False):
    """
    训练模型，评估并绘制结果。
    """
    print(f"\n--- {model_name} ---")
    mlp_model.fit(X_train, y_train)

    y_train_pred = mlp_model.predict(X_train)
    y_test_pred = mlp_model.predict(X_test)

    print(f"  训练集 R^2 Score: {r2_score(y_train, y_train_pred):.4f}")
    print(f"  测试集 R^2 Score: {r2_score(y_test, y_test_pred):.4f}")
    print(f"  测试集 MAE: {mean_absolute_error(y_test, y_test_pred):.4f}")
    print(f"  测试集 MSE: {mean_squared_error(y_test, y_test_pred):.4f}")

    if plot_loss and hasattr(mlp_model, 'loss_curve_'):
        plt.figure(figsize=(10, 6))
        plt.plot(mlp_model.loss_curve_, label=f'{model_name} Loss Curve')
        plt.title(f'{model_name} Training Loss Curve')
        plt.xlabel('Epochs')
        plt.ylabel('Loss (MSE)')
        plt.legend()
        plt.grid(True)
        plt.show()

    if plot_predictions:
        plt.figure(figsize=(10, 6))
        plt.scatter(y_test, y_test_pred, alpha=0.7, label='Predicted vs Actual')
        plt.plot([min(y_test.min(), y_test_pred.min()), max(y_test.max(), y_test_pred.max())], 
                 [min(y_test.min(), y_test_pred.min()), max(y_test.max(), y_test_pred.max())], 
                 '--', color='red', label='Ideal Fit')
        plt.xlabel('Actual PM2.5')
        plt.ylabel('Predicted PM2.5')
        plt.title(f'{model_name}: Actual vs. Predicted PM2.5 on Test Set')
        plt.legend()
        plt.grid(True)
        plt.show()
    return mlp_model


def main_regression_task():
    """
    主函数，执行PM2.5回归预测任务。
    """
    # 创建一个隐藏的Tkinter根窗口
    root = tk.Tk()
    root.withdraw() # 隐藏主窗口

    # 打开文件选择对话框
    csv_file_path = filedialog.askopenfilename(
        title="请选择空气质量监测数据CSV文件",
        filetypes=(("CSV 文件", "*.csv"), ("所有文件", "*.*"))
    )

    if not csv_file_path:
        print("用户未选择文件。程序将退出。")
        return

    print(f"用户选择的文件: {csv_file_path}")

    X, y, feature_names, scaler_X = load_and_preprocess_data(csv_file_path)

    if X is None or y is None:
        print("数据加载或预处理失败，程序终止。")
        return

    if len(y) < 20: # 检查是否有足够的数据进行划分和训练
        print(f"警告：数据量过小 (只有 {len(y)} 条有效记录)，模型训练效果可能不佳。")
        if len(y) < 2: # 无法划分
             print("错误：有效数据不足2条，无法进行训练测试划分。")
             return


    # 划分数据集
    test_split_size = 0.25
    if len(y) * test_split_size < 1: 
        test_split_size = 1 / len(y) if len(y) > 1 else 0 
        if test_split_size == 0 and len(y) == 1: 
            print("错误：只有一条有效数据，无法划分训练集和测试集。")
            return
        print(f"数据量较小，测试集比例调整为 {test_split_size:.2f}")


    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_split_size, random_state=42)

    common_params = {
        'solver': 'adam',
        'max_iter': 1000, 
        'random_state': 1,
        'learning_rate_init': 0.001,
        'early_stopping': True, 
        'n_iter_no_change': 20, 
        'tol': 1e-5 
    }

    # 模型1: ReLU, 100 hidden nodes
    mlp_relu_100 = MLPRegressor(
        hidden_layer_sizes=(100,), 
        activation='relu', 
        **common_params
    )
    model1 = train_evaluate_plot_model(X_train, X_test, y_train, y_test, 
                                       "Model 1 (ReLU, 100 nodes)", mlp_relu_100, plot_predictions=True)

    # 模型2: Logistic, 100 hidden nodes
    mlp_logistic_100 = MLPRegressor(
        hidden_layer_sizes=(100,),
        activation='logistic',
        **common_params
    )
    model2 = train_evaluate_plot_model(X_train, X_test, y_train, y_test,
                                       "Model 2 (Logistic, 100 nodes)", mlp_logistic_100)

    # 模型3: ReLU, hidden_layer_sizes=(100, 50)
    mlp_relu_100_50 = MLPRegressor(
        hidden_layer_sizes=(100, 50),
        activation='relu',
        **common_params
    )
    model3 = train_evaluate_plot_model(X_train, X_test, y_train, y_test,
                                       "Model 3 (ReLU, 100, 50 nodes)", mlp_relu_100_50)
    
    # 综合绘制两种激活函数下的损失函数变化曲线 (针对100个隐藏节点)
    if model1 and model2 and hasattr(model1, 'loss_curve_') and hasattr(model2, 'loss_curve_'):
        plt.figure(figsize=(12, 7))
        plt.plot(model1.loss_curve_, label='ReLU (100 nodes) Loss')
        plt.plot(model2.loss_curve_, label='Logistic (100 nodes) Loss')
        plt.title('Comparison of Training Loss Curves (100 Hidden Nodes)')
        plt.xlabel('Epochs')
        plt.ylabel('Loss (MSE)')
        plt.legend()
        plt.grid(True)
        plt.show()

    print("\n任务1完成。脚本已尝试训练和评估三个神经网络模型，并绘制了相关图表。")

if __name__ == "__main__":
    main_regression_task()