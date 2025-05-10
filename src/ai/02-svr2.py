import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
import warnings

# 忽略一些sklearn的FutureWarning，如果出现的话
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning) # Matplotlib 可能的 UserWarning

# --- 0. Matplotlib 中文显示配置 ---
try:
    #plt.rcParams['font.sans-serif'] = ['SimHei']  # 指定默认字体：黑体
    plt.rcParams['font.sans-serif'] = ['STHeiti']  # Specify default font: PingFang SC
    plt.rcParams['axes.unicode_minus'] = False  # 解决保存图像是负号'-'显示为方块的问题
except Exception as e:
    print(f"设置中文字体失败: {e}. 图表中的中文可能无法正确显示。")
    print("请确保已安装 SimHei 字体，或修改为系统中存在的其他中文字体（如 Microsoft YaHei, WenQuanYi Micro Hei 等）。")

# --- 1. 数据加载和初步处理函数 ---
def load_and_preprocess_data():
    """
    使用Tkinter让用户选择CSV文件，加载数据并进行预处理。
    预处理包括：
    - 处理'horsepower'列的'?'为NaN，转换为数值型，并用均值填充NaN。
    - 删除'MPG', 'weight', 'horsepower'中任何包含NaN的行。
    返回:
        X (pd.DataFrame): 特征 ('weight', 'horsepower')
        y (pd.Series): 目标变量 ('MPG')
        df_cleaned (pd.DataFrame): 清理后的完整DataFrame，用于选择x0
    """
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title="请选择 汽车MPG.CSV 文件",
        filetypes=[("CSV files", "*.csv")]
    )
    if not file_path:
        print("未选择文件，程序退出。")
        return None, None, None

    try:
        df = pd.read_csv(file_path)
        print(f"成功从 '{file_path}' 加载数据。")
    except Exception as e:
        print(f"读取文件 '{file_path}' 失败: {e}")
        return None, None, None

    if 'horsepower' in df.columns:
        df['horsepower'] = pd.to_numeric(df['horsepower'], errors='coerce')
        if df['horsepower'].isnull().sum() > 0:
            print(f"列 'horsepower' 中发现 {df['horsepower'].isnull().sum()} 个非数值/缺失值，将用该列均值填充。")
            horsepower_mean = df['horsepower'].mean()
            df['horsepower'].fillna(horsepower_mean, inplace=True)
            print(f"使用 'horsepower' 均值 {horsepower_mean:.2f} 填充缺失值。")
    else:
        print("错误：数据集中未找到 'horsepower' 列。程序退出。")
        return None, None, None

    required_cols = ['MPG', 'weight', 'horsepower']
    for col in required_cols:
        if col not in df.columns:
            print(f"错误：数据集中缺少必需列 '{col}'。程序退出。")
            return None, None, None

    initial_rows = len(df)
    df.dropna(subset=required_cols, inplace=True)
    rows_dropped = initial_rows - len(df)
    if rows_dropped > 0:
        print(f"因 'MPG', 'weight', 'horsepower' 列存在缺失值，已删除 {rows_dropped} 行。")

    if df.empty:
        print("错误：预处理后数据为空。请检查CSV文件内容。程序退出。")
        return None, None, None

    X = df[['weight', 'horsepower']].copy()
    y = df['MPG'].copy()
    
    print(f"数据加载和预处理完成。最终数据集包含 {len(df)} 条记录。")
    return X, y, df

# --- 2. 交叉验证、模型训练与绘图函数 ---
def train_eval_and_plot_cv(X_full, y_full, model_template, model_name, n_splits=2, svr_extra_params_to_plot=None):
    """
    执行n折交叉验证，训练模型，并绘制散点图。
    对于SVR，可以额外绘制不同超参数下的回归线。
    参数:
        X_full (pd.DataFrame): 全部特征数据
        y_full (pd.Series): 全部目标数据
        model_template: sklearn兼容的模型实例 (如LinearRegression(), Pipeline(SVR()))
        model_name (str): 模型名称，用于图表标题和打印输出
        n_splits (int): 交叉验证的折数
        svr_extra_params_to_plot (list of dicts, optional): 
            仅当 model_name 包含 "SVR" 时使用。
            每个字典包含 {'C': float, 'epsilon': float, 'label': str, 'color': str}
            用于绘制额外的SVR回归线。
    返回:
        all_actuals (np.array): 所有折叠中测试集的实际MPG值 (基于 model_template)
        all_predictions (np.array): 所有折叠中测试集的预测MPG值 (基于 model_template)
    """
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    
    fig_weight, axs_weight = plt.subplots(n_splits, 1, figsize=(14, 6 * n_splits), squeeze=False) # 增加宽度
    fig_hp, axs_hp = plt.subplots(n_splits, 1, figsize=(14, 6 * n_splits), squeeze=False) # 增加宽度

    fig_weight.suptitle(f'{model_name} - MPG vs 自重 (Weight) - {n_splits}折交叉验证', fontsize=16)
    fig_hp.suptitle(f'{model_name} - MPG vs 马力 (Horsepower) - {n_splits}折交叉验证', fontsize=16)

    fold_num = 0
    all_actuals_list = []
    all_predictions_list = [] # 这些是基于 model_template 的预测

    for train_index, test_index in kf.split(X_full):
        X_train_fold, X_test_fold = X_full.iloc[train_index], X_full.iloc[test_index]
        y_train_fold, y_test_fold = y_full.iloc[train_index], y_full.iloc[test_index]

        # 训练主模型 (model_template) 并获取其预测
        current_model = model_template 
        current_model.fit(X_train_fold, y_train_fold)
        y_pred_fold = current_model.predict(X_test_fold) # 这是主模型的预测

        all_actuals_list.extend(y_test_fold)
        all_predictions_list.extend(y_pred_fold)

        # --- 绘制 MPG vs Weight ---
        ax_w = axs_weight[fold_num, 0]
        ax_w.scatter(X_train_fold['weight'], y_train_fold, color='lightblue', label='训练集点', alpha=0.5, s=30)
        ax_w.scatter(X_test_fold['weight'], y_test_fold, color='green', label='测试集实际值', alpha=0.8, s=50, marker='o')
        ax_w.scatter(X_test_fold['weight'], y_pred_fold, color='red', label=f'{model_name} 测试集预测值', alpha=0.8, s=50, marker='x')

        weight_range = np.linspace(X_full['weight'].min(), X_full['weight'].max(), 100)
        hp_mean_train = X_train_fold['horsepower'].mean()
        X_plot_w = pd.DataFrame({'weight': weight_range, 'horsepower': hp_mean_train})
        
        # 绘制主模型的回归线
        y_plot_pred_w_main = current_model.predict(X_plot_w)
        ax_w.plot(weight_range, y_plot_pred_w_main, color='darkorange', linestyle='--', linewidth=2.5, label=f'{model_name} 回归线 (HP={hp_mean_train:.0f})')
        
        # 如果是SVR且有额外参数，绘制它们的回归线
        if "SVR" in model_name and svr_extra_params_to_plot:
            for params in svr_extra_params_to_plot:
                svr_variant_pipeline = Pipeline([
                    ('scaler', StandardScaler()),
                    ('svr', SVR(C=params['C'], epsilon=params['epsilon'], kernel='rbf'))
                ])
                svr_variant_pipeline.fit(X_train_fold, y_train_fold)
                y_plot_pred_w_variant = svr_variant_pipeline.predict(X_plot_w)
                ax_w.plot(weight_range, y_plot_pred_w_variant, color=params['color'], linestyle=':', linewidth=2, label=f'{params["label"]} (HP={hp_mean_train:.0f})')

        ax_w.set_title(f'第 {fold_num + 1} 折', fontsize=12)
        ax_w.set_xlabel('自重 (Weight)', fontsize=10)
        ax_w.set_ylabel('MPG', fontsize=10)
        ax_w.legend(fontsize=8, loc='upper right')
        ax_w.grid(True, linestyle=':', alpha=0.7)

        # --- 绘制 MPG vs Horsepower ---
        ax_h = axs_hp[fold_num, 0]
        ax_h.scatter(X_train_fold['horsepower'], y_train_fold, color='lightblue', label='训练集点', alpha=0.5, s=30)
        ax_h.scatter(X_test_fold['horsepower'], y_test_fold, color='green', label='测试集实际值', alpha=0.8, s=50, marker='o')
        ax_h.scatter(X_test_fold['horsepower'], y_pred_fold, color='red', label=f'{model_name} 测试集预测值', alpha=0.8, s=50, marker='x')

        hp_range = np.linspace(X_full['horsepower'].min(), X_full['horsepower'].max(), 100)
        weight_mean_train = X_train_fold['weight'].mean()
        X_plot_h = pd.DataFrame({'weight': weight_mean_train, 'horsepower': hp_range})

        # 绘制主模型的回归线
        y_plot_pred_h_main = current_model.predict(X_plot_h)
        ax_h.plot(hp_range, y_plot_pred_h_main, color='purple', linestyle='--', linewidth=2.5, label=f'{model_name} 回归线 (W={weight_mean_train:.0f})')

        # 如果是SVR且有额外参数，绘制它们的回归线
        if "SVR" in model_name and svr_extra_params_to_plot:
            for params in svr_extra_params_to_plot:
                svr_variant_pipeline = Pipeline([
                    ('scaler', StandardScaler()),
                    ('svr', SVR(C=params['C'], epsilon=params['epsilon'], kernel='rbf'))
                ])
                svr_variant_pipeline.fit(X_train_fold, y_train_fold)
                y_plot_pred_h_variant = svr_variant_pipeline.predict(X_plot_h)
                ax_h.plot(hp_range, y_plot_pred_h_variant, color=params['color'], linestyle=':', linewidth=2, label=f'{params["label"]} (W={weight_mean_train:.0f})')
        
        ax_h.set_title(f'第 {fold_num + 1} 折', fontsize=12)
        ax_h.set_xlabel('马力 (Horsepower)', fontsize=10)
        ax_h.set_ylabel('MPG', fontsize=10)
        ax_h.legend(fontsize=8, loc='upper right')
        ax_h.grid(True, linestyle=':', alpha=0.7)
        
        # 以下指标基于主模型 (model_template)
        pred_mean_fold = np.mean(y_pred_fold)
        pred_var_fold = np.var(y_pred_fold)
        error_fold = y_test_fold.values - y_pred_fold
        error_mean_fold = np.mean(error_fold)
        error_var_fold = np.var(error_fold)
        r2_fold = r2_score(y_test_fold, y_pred_fold)
        mse_fold = mean_squared_error(y_test_fold, y_pred_fold)

        print(f"  {model_name} (主模型) - 第 {fold_num + 1} 折:")
        print(f"    测试集预测MPG: 均值={pred_mean_fold:.2f}, 方差={pred_var_fold:.2f}")
        print(f"    测试集预测误差: 均值={error_mean_fold:.2f}, 方差={error_var_fold:.2f}")
        print(f"    测试集R²得分: {r2_fold:.3f}, MSE: {mse_fold:.2f}")

        fold_num += 1

    fig_weight.tight_layout(rect=[0, 0, 1, 0.96]) 
    fig_hp.tight_layout(rect=[0, 0, 1, 0.96])

    return np.array(all_actuals_list), np.array(all_predictions_list)


# --- 3. 特定数据点 x0 的预测均值和方差函数 ---
def predict_for_x0(X_full, y_full, x0_features_df, model_template, model_name, n_splits=2):
    """
    对指定数据点x0，使用交叉验证中每折训练的主模型进行预测，并计算预测的均值和方差。
    """
    predictions_on_x0 = []
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42) 

    for train_index, _ in kf.split(X_full): 
        X_train_fold, y_train_fold = X_full.iloc[train_index], y_full.iloc[train_index]
        
        # 确保每次都用模板创建一个新的模型实例或正确地重新训练
        # 对于Pipeline，直接fit会重新训练内部的转换器和估计器
        # 对于简单模型如LinearRegression，fit会重置并训练
        current_model_instance = model_template # 使用传递进来的主模型模板
        current_model_instance.fit(X_train_fold, y_train_fold)
        
        pred_x0 = current_model_instance.predict(x0_features_df)[0]
        predictions_on_x0.append(pred_x0)

    mean_pred_x0 = np.mean(predictions_on_x0)
    var_pred_x0 = np.var(predictions_on_x0)

    print(f"\n--- 对 x0 ({x0_features_df.iloc[0].to_dict()}) 使用 {model_name} (主模型) 的预测 ({n_splits}折模型) ---")
    if not predictions_on_x0:
        print("  未能生成任何预测值。")
        return np.nan, np.nan
        
    print(f"  各折模型对 x0 的预测值: {[f'{p:.2f}' for p in predictions_on_x0]}")
    print(f"  对 x0 的预测均值: {mean_pred_x0:.2f} MPG")
    print(f"  对 x0 的预测方差: {var_pred_x0:.2f}")
    return mean_pred_x0, var_pred_x0

# --- 4. 主程序 ---
def main():
    X, y, df_cleaned = load_and_preprocess_data()
    if X is None or y is None or df_cleaned is None:
        return 

    # 定义主模型
    lr_model_template = LinearRegression()
    # 定义SVR的主模型 (用于性能比较和x0预测)
    default_svr_C = 1.0
    default_svr_epsilon = 0.1
    svr_model_template = Pipeline([
        ('scaler', StandardScaler()),  
        ('svr', SVR(kernel='rbf', C=default_svr_C, epsilon=default_svr_epsilon))      
    ])
    
    n_cv_splits = 2 

    print(f"\n=== 开始执行 {n_cv_splits}折交叉验证：线性回归 ===")
    lr_actuals, lr_predictions = train_eval_and_plot_cv(
        X, y, lr_model_template, "线性回归 (LR)", n_splits=n_cv_splits
    )
    overall_lr_mse = mean_squared_error(lr_actuals, lr_predictions)
    overall_lr_r2 = r2_score(lr_actuals, lr_predictions)
    print(f"线性回归 ({n_cv_splits}-fold CV 总体性能):")
    print(f"  MSE: {overall_lr_mse:.2f}")
    print(f"  R² score: {overall_lr_r2:.3f}")

    # 定义用于额外绘制的SVR参数组合
    svr_extra_params_for_plotting = [
        {'C': 0.1, 'epsilon': 0.1, 'label': 'SVR (C=0.1, ε=0.1)', 'color': 'cyan'},
        {'C': 10, 'epsilon': 0.1, 'label': 'SVR (C=10, ε=0.1)', 'color': 'magenta'},
        {'C': default_svr_C, 'epsilon': 0.5, 'label': f'SVR (C={default_svr_C}, ε=0.5)', 'color': 'lime'},
        # 可以添加更多组合
    ]
    # 确保默认SVR的标签是唯一的，或者主回归线标签已足够
    main_svr_label_part = f"SVR (主: C={default_svr_C}, ε={default_svr_epsilon})"

    print(f"\n=== 开始执行 {n_cv_splits}折交叉验证：支持向量回归 ===")
    print(f"SVR主模型参数: C={default_svr_C}, epsilon={default_svr_epsilon}")
    print("将额外绘制以下SVR参数组合的回归线:")
    for p in svr_extra_params_for_plotting:
        print(f"  C={p['C']}, epsilon={p['epsilon']}")

    svr_actuals, svr_predictions = train_eval_and_plot_cv(
        X, y, svr_model_template, main_svr_label_part, 
        n_splits=n_cv_splits, 
        svr_extra_params_to_plot=svr_extra_params_for_plotting
    )
    overall_svr_mse = mean_squared_error(svr_actuals, svr_predictions) # 基于主SVR模型
    overall_svr_r2 = r2_score(svr_actuals, svr_predictions) # 基于主SVR模型
    print(f"支持向量回归 (主模型: C={default_svr_C}, ε={default_svr_epsilon}) ({n_cv_splits}-fold CV 总体性能):")
    print(f"  MSE: {overall_svr_mse:.2f}")
    print(f"  R² score: {overall_svr_r2:.3f}")

    print("\n=== 模型比较总结 (基于主模型整体交叉验证结果) ===")
    print(f"线性回归 (LR):     MSE = {overall_lr_mse:.2f}, R² = {overall_lr_r2:.3f}")
    print(f"支持向量回归 (SVR 主模型 C={default_svr_C}, ε={default_svr_epsilon}): MSE = {overall_svr_mse:.2f}, R² = {overall_svr_r2:.3f}")
    if overall_lr_mse < overall_svr_mse:
        print("在此数据集和评估下，线性回归的MSE较低，表现相对较好。")
    elif overall_svr_mse < overall_lr_mse:
        print("在此数据集和评估下，支持向量回归(主模型)的MSE较低，表现相对较好。")
    else:
        print("在此数据集和评估下，两种模型(主模型)的MSE表现相似。")

    # --- 对指定数据点 x0 进行预测 (使用主模型) ---
    if len(X) > 0:
        idx_pos_x0 = np.random.randint(0, len(X)) 
        print(f"\n随机选择数据行（基于处理后的X/y数据集的位置索引 {idx_pos_x0}）作为 x0。")
        
        x0_features_df = X.iloc[[idx_pos_x0]]
        x0_actual_mpg = y.iloc[idx_pos_x0]    
        
        original_index_label = X.index[idx_pos_x0] 
        print(f"选定的 x0 (原始DataFrame索引标签: {original_index_label}): 特征 {x0_features_df.iloc[0].to_dict()}")
        print(f"x0 的实际 MPG: {x0_actual_mpg:.2f}")

        predict_for_x0(X, y, x0_features_df, lr_model_template, "线性回归 (LR)", n_splits=n_cv_splits)
        predict_for_x0(X, y, x0_features_df, svr_model_template, f"支持向量回归 (主模型 C={default_svr_C}, ε={default_svr_epsilon})", n_splits=n_cv_splits)
    else:
        print("数据为空，无法选择 x0 进行预测。")

    print("\n所有图表已生成。关闭图表窗口后程序结束。")
    plt.show()

if __name__ == '__main__':
    main()