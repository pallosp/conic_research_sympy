import contextlib
import io
import sys
from pathlib import Path

import yaml
from pydoc_markdown import PydocMarkdown
from pydoc_markdown.interfaces import Context

CONFIG_FILE = Path("pydoc-markdown.yml")
API_FILE = Path("docs/api.md")
LIB_DIR = Path("src/lib")


def main() -> None:
    """Ensures that docs/api.md is up to date with pydoc-markdown."""
    config = yaml.safe_load(CONFIG_FILE.read_text(encoding="utf-8"))

    # The config lists the modules explicitly to control their order in the docs.
    listed = config["loaders"][0]["modules"]
    existing = [p.stem for p in LIB_DIR.glob("*.py") if p.stem != "__init__"]
    if sorted(listed) != sorted(existing):
        sys.exit(f"{CONFIG_FILE} must list each module of {LIB_DIR} exactly once.")

    # Render to stdout instead of the configured file.
    del config["renderer"]["filename"]
    session = PydocMarkdown()
    session.load_config(config)
    # Needed for the source links. The directory must be absolute.
    session.init(Context(directory=str(CONFIG_FILE.parent.absolute())))
    modules = session.load_modules()
    session.process(modules)
    generated = io.StringIO()
    with contextlib.redirect_stdout(generated):
        session.render(modules)

    if generated.getvalue() != API_FILE.read_text(encoding="utf-8"):
        sys.exit(f"{API_FILE} is out of date. Run 'pydoc-markdown' to regenerate.")


if __name__ == "__main__":
    main()
