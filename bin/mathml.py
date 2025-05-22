import re
import sys
import pyperclip
from markdown import markdown
import latex2mathml.converter

def convert_math_to_mathml(text):
    # Convert inline math $...$
    inline_math_pattern = re.compile(r'\$(.*?)\$')
    text = inline_math_pattern.sub(lambda m: convert_to_mathml(m.group(1)), text)

    # Convert block math $$...$$
    block_math_pattern = re.compile(r'\$\$(.*?)\$\$', re.DOTALL)
    text = block_math_pattern.sub(lambda m: convert_to_mathml(m.group(1)), text)

    return text

def convert_to_mathml(latex_expr):
    try:
        # Convert LaTeX to MathML using latex2mathml
        mathml = latex2mathml.converter.convert(latex_expr)
        return mathml
    except Exception as e:
        print(f"Error converting LaTeX to MathML: {e}")
        return latex_expr

def process_markdown_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        content = file.read()

    # Convert Markdown to HTML
    html_content = markdown(content)

    # Convert math expressions to MathML
    converted_content = convert_math_to_mathml(html_content)

    # Copy converted content to clipboard
    pyperclip.copy(converted_content)
    print("Converted content copied to clipboard.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <markdown_file>")
        sys.exit(1)

    markdown_file = sys.argv[1]
    process_markdown_file(markdown_file)