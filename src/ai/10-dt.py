import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor, export_text, plot_tree
from sklearn.metrics import accuracy_score, mean_squared_error
import matplotlib.pyplot as plt # 统一导入pyplot
from mpl_toolkits.mplot3d import Axes3D # 用于3D绘图
from io import StringIO
import tkinter as tk
from tkinter import filedialog

# --- Matplotlib 中文显示设置 ---
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'STHeiti']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False

# 提供的示例CSV数据，作为备用
EXAMPLE_CSV_DATA = """期,AQI,质量等级,PM2.5,PM10,SO2,CO,NO2,O3
1/1/14,81,良,45,111,28,1.5,62,52
1/2/14,145,轻度污染,111,168,69,3.4,93,14
1/3/14,74,良,47,98,29,1.3,52,56
1/4/14,149,轻度污染,114,147,40,2.8,75,14
1/5/14,119,轻度污染,91,117,36,2.3,67,44
1/6/14,182,中度污染,138,158,46,2.4,68,12
1/7/14,145,轻度污染,111,125,34,2,60,43
1/8/14,27,优,15,25,13,0.5,21,53
1/9/14,46,优,27,46,19,0.8,35,53
1/10/14,85,良,63,94,53,1.9,71,19
1/11/14,139,轻度污染,106,128,76,2.8,90,11
1/12/14,47,优,27,47,27,0.7,39,59
1/13/14,109,轻度污染,82,107,67,2.3,78,20
1/14/14,108,轻度污染,82,108,68,2.4,74,24
1/15/14,125,轻度污染,95,117,66,2.7,84,18
1/16/14,402,严重污染,353,384,109,4.6,123,20
1/17/14,215,重度污染,165,194,67,2.4,76,30
1/18/14,99,良,74,81,42,1.5,61,36
1/19/14,151,中度污染,115,186,71,2.5,72,41
"""

def preprocess_data(df):
    """对DataFrame进行预处理：清理列名、转换数据类型、处理NaN"""
    if df is None or df.empty:
        return None
    df.columns = df.columns.str.strip()
    cols_to_numeric = ['PM2.5', 'PM10', 'SO2', 'CO', 'NO2', 'AQI', 'O3']
    existing_cols_to_numeric = [col for col in cols_to_numeric if col in df.columns]
    for col in existing_cols_to_numeric:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    if df[existing_cols_to_numeric].isnull().any().any():
        print("警告: 数值列转换后存在NaN值。为继续运行，将用0填充这些NaN值。")
        for col in existing_cols_to_numeric:
            if df[col].isnull().any():
                df[col].fillna(0, inplace=True)
    return df

def load_data_from_filepath(filepath):
    """从指定文件路径加载并预处理CSV数据"""
    try:
        df = pd.read_csv(filepath)
        return preprocess_data(df)
    except FileNotFoundError:
        print(f"错误: 文件未找到: {filepath}")
        return None
    except pd.errors.EmptyDataError:
        print(f"错误: 文件为空: {filepath}")
        return None
    except Exception as e:
        print(f"读取或预处理文件 {filepath} 时发生错误: {e}")
        return None

def load_data_from_string(csv_string):
    """从字符串加载并预处理CSV数据"""
    try:
        df = pd.read_csv(StringIO(csv_string))
        return preprocess_data(df)
    except Exception as e:
        print(f"从内置字符串数据加载或预处理时发生错误: {e}")
        return None

def classification_task(df_orig):
    """执行分类任务"""
    if df_orig is None or df_orig.empty:
        print("分类任务错误：输入数据为空。")
        return
    df = df_orig.copy()
    print("\n--- 开始执行分类任务 ---")
    
    features_classification = ['PM2.5', 'PM10', 'SO2', 'CO', 'NO2']
    target_classification = '质量等级'

    missing_cols = [col for col in features_classification + [target_classification] if col not in df.columns]
    if missing_cols:
        print(f"错误: 数据中缺少以下列: {', '.join(missing_cols)}。无法执行分类任务。")
        return

    X = df[features_classification]
    y = df[target_classification]

    if len(df) < 5:
        print("数据量过少（少于5行），无法有效进行训练集和测试集分割及模型评估。")
        if len(df) > 0 : 
            print("将尝试仅用全部数据训练模型并展示其结构。")
            clf_small_data = DecisionTreeClassifier(criterion='entropy', random_state=42)
            try:
                clf_small_data.fit(X,y)
                print("\n分类树 (文本化 - 基于全部少量数据):")
                tree_text_small = export_text(clf_small_data, feature_names=list(X.columns))
                print(tree_text_small)
                
                plt.figure(figsize=(15, 8))
                plot_tree(clf_small_data, feature_names=list(X.columns), class_names=sorted(list(y.unique())), filled=True, rounded=True, fontsize=8)
                plt.title("分类决策树 (基于全部少量数据)")
                print("少量数据下的分类树图形即将显示。您可以使用绘图窗口的保存按钮保存图像。")
                plt.show()
            except Exception as e:
                print(f"在少量数据上训练或绘图时发生错误: {e}")
        return

    X_train, X_test, y_train, y_test = pd.DataFrame(), pd.DataFrame(), pd.Series(dtype='object'), pd.Series(dtype='object')
    try:
        if y.nunique() > 1 and all(y.value_counts() >= 2):
             X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
        else:
            print("警告: 类别数量不足或某些类别样本过少，无法进行分层抽样。将使用普通随机抽样。")
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    except ValueError as e:
        print(f"警告: 划分数据集时发生错误 ({e})。将使用普通随机抽样。")
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    if X_train.empty or y_train.empty:
        print("错误：训练集为空，无法继续分类任务。")
        return

    clf = DecisionTreeClassifier(criterion='entropy', random_state=42) 
    clf.fit(X_train, y_train)

    y_train_pred = clf.predict(X_train)
    train_accuracy = accuracy_score(y_train, y_train_pred)
    print(f"训练精度: {train_accuracy:.4f}")

    if not X_test.empty:
        y_test_pred = clf.predict(X_test)
        test_accuracy = accuracy_score(y_test, y_test_pred)
        print(f"测试精度: {test_accuracy:.4f}")
    else:
        print("测试集为空，无法计算测试精度。")

    if train_accuracy < 0.80:
        print("提示: 训练精度低于要求的80%。")

    print("\n分类树 (文本化):")
    feature_names_list = list(X.columns)
    class_names_list = sorted(list(y.unique())) 
    tree_text = export_text(clf, feature_names=feature_names_list)
    print(tree_text)

    print("\n生成分类树图形...")
    plt.figure(figsize=(20, 12))
    plot_tree(clf,
              feature_names=feature_names_list,
              class_names=class_names_list,
              filled=True,
              rounded=True,
              fontsize=7) 
    plt.title("分类决策树 (criterion='entropy')")
    print("分类树图形即将显示。您可以使用绘图窗口的保存按钮保存图像。")
    plt.show()


def regression_task(df_orig):
    """执行回归任务"""
    if df_orig is None or df_orig.empty:
        print("回归任务错误：输入数据为空。")
        return
    df = df_orig.copy()
    print("\n--- 开始执行回归任务 ---")

    features_regression = ['CO', 'SO2']
    target_regression = 'PM2.5'

    missing_cols_reg = [col for col in features_regression + [target_regression] if col not in df.columns]
    if missing_cols_reg:
        print(f"错误: 数据中缺少以下列: {', '.join(missing_cols_reg)}。无法执行回归任务。")
        return
        
    X = df[features_regression]
    y = df[target_regression]

    if len(df) < 5:
        print("数据量过少（少于5行），回归任务的误差曲线和最优深度可能不具代表性。")
        return

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    if X_train.empty or y_train.empty:
        print("错误：训练集为空，无法继续回归任务。")
        return

    depths = range(2, 16)
    train_errors = []
    test_errors = []

    for depth in depths:
        reg = DecisionTreeRegressor(max_depth=depth, random_state=42)
        reg.fit(X_train, y_train)
        y_train_pred = reg.predict(X_train)
        train_errors.append(mean_squared_error(y_train, y_train_pred))
        if not X_test.empty:
            y_test_pred = reg.predict(X_test)
            test_errors.append(mean_squared_error(y_test, y_test_pred))
        else:
            test_errors.append(float('nan'))

    plt.figure(figsize=(12, 7))
    plt.plot(depths, train_errors, marker='o', linestyle='-', label='训练误差 (MSE)')
    if not X_test.empty and any(not np.isnan(err) for err in test_errors):
        plt.plot(depths, test_errors, marker='s', linestyle='--', label='测试误差 (MSE)')
    plt.xlabel('树深度 (max_depth)')
    plt.ylabel('均方误差 (MSE)')
    plt.title('回归树深度与误差关系图')
    plt.legend()
    plt.grid(True)
    plt.xticks(list(depths))
    print("树深度与误差折线图即将显示。您可以使用绘图窗口的保存按钮保存图像。")
    plt.show()

    optimal_depth = -1
    min_test_error_val = float('inf')

    if not X_test.empty and any(not np.isnan(err) for err in test_errors):
        valid_test_errors = [(d, err) for d, err in zip(depths, test_errors) if not np.isnan(err)]
        if valid_test_errors:
            min_test_error_val = min(err for d, err in valid_test_errors)
            optimal_depth_candidates = [d for d, err in valid_test_errors if np.isclose(err, min_test_error_val)]
            if optimal_depth_candidates:
                optimal_depth = min(optimal_depth_candidates)
                print(f"\n基于测试误差，选定的最优树深度为: {optimal_depth} (最小测试MSE: {min_test_error_val:.4f})")
    
    if optimal_depth == -1 : 
        print("未能通过测试误差有效确定最优深度。")
        if train_errors:
            min_train_error_val = min(train_errors)
            optimal_depth_train_candidates = [depths[i] for i, err in enumerate(train_errors) if np.isclose(err, min_train_error_val)]
            if optimal_depth_train_candidates:
                 optimal_depth = min(optimal_depth_train_candidates)
                 print(f"备选方案：基于训练误差最小选择的深度为: {optimal_depth} (训练MSE: {min_train_error_val:.4f})")
                 print("注意：此深度可能导致过拟合。")
            else:
                 optimal_depth = 3 
                 print(f"无法从训练误差中确定备选深度。将使用默认深度 {optimal_depth}。")
        else:
            optimal_depth = 3
            print(f"无有效的训练误差数据。将使用默认深度 {optimal_depth}。")
    
    if optimal_depth == -1:
        optimal_depth = 3
        print(f"无法确定最优深度，将使用默认深度 {optimal_depth} 进行后续绘图。")

    print(f"最终选定的最优树深度 (用于绘图): {optimal_depth}")

    optimal_reg = DecisionTreeRegressor(max_depth=optimal_depth, random_state=42)
    optimal_reg.fit(X_train, y_train)

    # --- 新增：绘制最优回归树的图形 ---
    print("\n生成最优回归树图形...")
    plt.figure(figsize=(18, 10)) # 可以根据树的复杂度调整大小
    plot_tree(optimal_reg,
              feature_names=features_regression, # 回归任务的特征名
              filled=True,
              rounded=True,
              fontsize=7,
              precision=2) # 显示回归值的精度
    plt.title(f"最优回归决策树 (深度: {optimal_depth})")
    print(f"最优回归树 (深度 {optimal_depth}) 图形即将显示。您可以使用绘图窗口的保存按钮保存图像。")
    plt.show()
    # --- 新增结束 ---

    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')

    co_min_overall, co_max_overall = X[features_regression[0]].min(), X[features_regression[0]].max()
    so2_min_overall, so2_max_overall = X[features_regression[1]].min(), X[features_regression[1]].max()
    
    co_extra = (co_max_overall - co_min_overall) * 0.1 if (co_max_overall - co_min_overall) > 0 else 0.5
    so2_extra = (so2_max_overall - so2_min_overall) * 0.1 if (so2_max_overall - so2_min_overall) > 0 else 0.5

    co_mesh_vals = np.linspace(co_min_overall - co_extra, co_max_overall + co_extra, 30)
    so2_mesh_vals = np.linspace(so2_min_overall - so2_extra, so2_max_overall + so2_extra, 30)
    
    CO_mesh, SO2_mesh = np.meshgrid(co_mesh_vals, so2_mesh_vals)
    
    mesh_features_data = np.c_[CO_mesh.ravel(), SO2_mesh.ravel()]
    mesh_features_df = pd.DataFrame(mesh_features_data, columns=features_regression)
    PM25_pred_mesh = optimal_reg.predict(mesh_features_df).reshape(CO_mesh.shape)

    ax.plot_surface(CO_mesh, SO2_mesh, PM25_pred_mesh, alpha=0.5, cmap='viridis', rstride=1, cstride=1, edgecolor='none')

    if not X_test.empty and not y_test.empty:
        y_pred_test_optimal = optimal_reg.predict(X_test)
        colors = ['blue' if pred < actual else 'grey' 
                  for pred, actual in zip(y_pred_test_optimal, y_test)]
        
        ax.scatter(X_test[features_regression[0]], X_test[features_regression[1]], y_test, 
                   c=colors, marker='o', s=60, depthshade=True, label='实际PM2.5 (测试集)')
        
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', label='预测值 < 实际值', markersize=8, markerfacecolor='blue'),
            Line2D([0], [0], marker='o', color='w', label='预测值 >= 实际值', markersize=8, markerfacecolor='grey')
        ]
        ax.legend(handles=legend_elements, title="测试点颜色", loc="upper left")
    else:
        ax.legend(loc="upper left")

    ax.set_xlabel(f'{features_regression[0]} 浓度')
    ax.set_ylabel(f'{features_regression[1]} 浓度')
    ax.set_zlabel(f'{target_regression} 浓度')
    ax.set_title(f'最优树深度 ({optimal_depth})下的 PM2.5 回归平面')
    print("回归平面图即将显示。您可以使用绘图窗口的保存按钮保存图像。")
    plt.show()


def main():
    """主函数，执行实验任务"""
    df = None
    root = tk.Tk()
    root.withdraw() 
    filepath = filedialog.askopenfilename(
        title="请选择一个CSV空气质量数据文件",
        filetypes=(("CSV 文件", "*.csv"), ("所有文件", "*.*"))
    )
    
    if filepath:
        print(f"尝试从文件加载数据: {filepath}")
        df = load_data_from_filepath(filepath)
    else:
        print("未选择文件。")

    if df is None or df.empty:
        print("无法从文件加载数据或文件未选择。将使用内置的示例数据。")
        df = load_data_from_string(EXAMPLE_CSV_DATA)

    if df is None or df.empty:
        print("错误：无法加载任何数据（文件或示例）。程序将退出。")
        return

    classification_task(df)
    regression_task(df)

    print("\n--- 实验任务完成 ---")
    print("所有图形已按顺序显示。您可以使用各个绘图窗口的保存功能来保存图像。")

if __name__ == "__main__":
    main()