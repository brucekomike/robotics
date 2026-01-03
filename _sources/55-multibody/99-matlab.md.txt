# matlab/octave 建模
## 代码习惯
代码开头加入
```matlab
% 清空变量
clear all
% 清空命令窗口
clc
```

## octave 拓展
octave需要额外的配置以使用解析解
```matlab
usepackage('symbolic');
```