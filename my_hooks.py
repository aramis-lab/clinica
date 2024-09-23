def on_startup(markdown, *, page, config, files):
    return markdown.replace("__DOCS_DIR__", "startup")


# def on_serve(markdown, *, page, config, files):
#     return markdown.replace("__DOCS_DIR__", "serve")
#
# def on_page_markdown(markdown, *, page, config, files):
#     return markdown.replace("__DOCS_DIR__", "page")
#
# def on_pre_build(markdown, *, page, config, files):
#     return markdown.replace("__DOCS_DIR__", "pre_build")
