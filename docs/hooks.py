def on_page_markdown(markdown, *, page, config, files):
    return markdown.replace("__DOCS_DIR__", config.docs_dir)
