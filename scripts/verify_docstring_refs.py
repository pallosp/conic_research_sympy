import ast
import re
import sys
from collections.abc import Iterator
from pathlib import Path

# Regex for docstring references: [text](#anchor)
REF_PATTERN = re.compile(r"\[[^\]]+\]\(#([a-zA-Z0-9_\.]+)\)")

# AST nodes that have a 'name' attribute and can be referenced
NAMED_TYPES = (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
# AST nodes that can have docstrings
DOC_TYPES = (ast.Module, *NAMED_TYPES)


def extract_anchors_from_tree(module: str, tree: ast.Module) -> Iterator[str]:
    """Yields all available documentation anchors from a single parsed tree."""
    yield module
    for node in tree.body:
        if isinstance(node, NAMED_TYPES):
            yield f"{module}.{node.name}"
            if isinstance(node, ast.ClassDef):
                for sub in node.body:
                    if isinstance(sub, NAMED_TYPES):
                        yield f"{module}.{node.name}.{sub.name}"
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name):
                    yield f"{module}.{target.id}"


def collect_anchors(trees: dict[Path, ast.Module]) -> Iterator[str]:
    """Yields all available documentation anchors from the parsed trees."""
    for f, tree in trees.items():
        if f.stem != "__init__":
            yield from extract_anchors_from_tree(f.stem, tree)


def main() -> None:
    """Main entry point for docstring reference verification."""
    src_root = Path("src/lib")
    files = list(src_root.glob("**/*.py"))

    # Parse all files once. Fails immediately on SyntaxError or OSError.
    trees = {
        f: ast.parse(f.read_text(encoding="utf-8"), filename=str(f)) for f in files
    }

    anchors = set(collect_anchors(trees))
    errors = []

    for f, tree in trees.items():
        for node in ast.walk(tree):
            if isinstance(node, DOC_TYPES):
                doc = ast.get_docstring(node)
                if not doc:
                    continue
                for ref in REF_PATTERN.findall(doc):
                    if ref not in anchors:
                        name = getattr(node, "name", f"module {f.stem}")
                        errors.append(f"{f}: In '{name}': Broken reference '{ref}'")

    if errors:
        for e in sorted(errors):
            print(e, file=sys.stderr)
        sys.exit(1)

    print("All docstring references verified successfully.")


if __name__ == "__main__":
    main()
