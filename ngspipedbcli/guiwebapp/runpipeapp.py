from click_web import create_click_web_app

import runpipe_click

app = create_click_web_app(runpipe_click, runpipe_click.cli)