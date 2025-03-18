

# MyST 深度扩展

## 1. 核心扩展功能

### 内联分页标签
````markdown
```{tab-set}
```{tab-item} Python
print("分页内容 1")
```

```{tab-item} C++
cout << "分页内容 2" << endl;
```

```{tab-item} Java
System.out.println("分页内容 3");
```
```
````
**渲染效果：** 生成可切换的代码分页

### 文献引用系统
1. 创建 `refs.bib` 文件：
```bibtex
@book{knuth1984tex,
  title={The TeXbook},
  author={Knuth, Donald E.},
  year={1984},
  publisher={Addison-Wesley}
}
```

2. 文档中引用：
```markdown
根据 {cite}`knuth1984tex` 的理论...
```

3. 添加参考文献章节：
```{bibliography}
:filter: cited
```

## 2. 高级指令系统

### 可折叠区块
```markdown
```{toggle}
隐藏的附加说明...
```
```

### 侧边注释栏
```markdown
```{sidebar} 重要提示
这是侧边栏内容
```
```

### 进度指示器
```markdown
```{progress} 65
当前进度 65%
```
```

## 3. 交互式组件

### 下拉菜单
```markdown
```{dropdown} 点击展开
隐藏的详细内容...
```
```

### 选择题
```markdown
```{choice} question1
:correct: 2

1. 错误选项
2. 正确答案
```
```

## 4. 科学文档支持

### 定理环境
```markdown
```{theorem} 勾股定理
:label: my-theorem

直角三角形斜边平方等于两直角边平方和
```
```

### 化学公式
```markdown
```{chemistry}
CH_3COOH ⇌ H^+ + CH_3COO^-
```
```

## 5. 多媒体集成

### 视频嵌入
```markdown
```{video} assets/demo.mp4
:width: 600
:poster: poster.jpg
```
```

### PDF 嵌入
```markdown
```{pdf} docs/spec.pdf
:page: 3
:width: 100%
```
```

## 6. 高级引用系统

### 跨文档引用
```markdown
参见 {doc}`../intro` 中的说明
```

### 自动编号公式
```markdown
```{math}
:label: eq1
E = mc^2
```

引用公式 {eq}`eq1`
```

## 7. 自定义扩展

### 添加插件
在 `_config.yml` 中配置：
```yaml
myst_extensions:
  - colon_fence
  - substitutions
  - deflist
```

### 自定义指令
创建 `custom.py`：
```python
from myst_parser.parsers.docutils import Directive

class Version(Directive):
    def run(self):
        return [nodes.strong("当前版本：2.1.0")]
```

## 8. 扩展对比表

| 功能             | 标准 Markdown | MyST 扩展 |
|------------------|---------------|-----------|
| 分页切换         | ❌            | ✅        |
| 文献引用         | ❌            | ✅        |
| 交互式组件       | ❌            | ✅        |
| 定理环境         | ❌            | ✅        |
| 跨文档引用       | ❌            | ✅        |
| 自定义指令       | ❌            | ✅        |

```{seealso}
- [furo theme](https://pradyunsg.me/furo)
- [myst](https://myst.tools/)
- [sphinx](https://www.sphinx-doc.org/)
```

> 提示：启用扩展需在配置文件 `_config.yml` 中声明