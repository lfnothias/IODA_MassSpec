from IPython.display import FileLink

class DownloadFileLink(FileLink):
    html_link_str = "<a href='{link}' download={file_name}>{link_text}</a>"

    def __init__(self, path, file_name=None, link_text=None, *args, **kwargs):
        super(DownloadFileLink, self).__init__(path, *args, **kwargs)

        self.file_name = file_name or os.path.split(path)[1]
        self.link_text = link_text or self.file_name

    def _format_path(self):
        from html import escape
        fp = ''.join([self.url_prefix, escape(self.path)])
        return ''.join([self.result_html_prefix,
                        self.html_link_str.format(link=fp, file_name=self.file_name, link_text=self.link_text),
                        self.result_html_suffix])
    