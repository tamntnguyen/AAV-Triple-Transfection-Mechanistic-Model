"""
Standard Dash dashboard navigation bar
"""

import dash_mantine_components as dmc
from dash import html, dcc, Dash, Output, Input, State
from dash_iconify import DashIconify

# DEVELOPER MODIFIED CONSTANTS
# ============================
# Change these to fit your needs for the dashboard
BAR_TITLE = "AAV Replication and Production Modeling"
TAB_TITLE = "AAV ODE Model"
GIT_URL = "https://github.com/tamntnguyen/AAV-Triple-Transfection-Mechanistic-Model"

README_MARKDOWN_FILE = "assets/README.md"  # Location of README markdown

# Set up the app
app = Dash(__name__)
app.title = TAB_TITLE


def create_home_link(label):
    """Create url

    Create a URL with name `label`

    Args:
      label(str): Label
    """
    return dmc.Text(
        label,
        size="xl",
        color="gray",
    )


GIT = html.A(
    dmc.Tooltip(
        dmc.ThemeIcon(
            DashIconify(
                icon="radix-icons:github-logo",
                width=22,
            ),
            radius=30,
            size=36,
            variant="outline",
            color="gray",
        ),
        label="Source Code",
        position="bottom",
    ),
    href=GIT_URL,
)

MORE_INFO = dmc.ActionIcon(
    dmc.Tooltip(
        dmc.ThemeIcon(
            DashIconify(
                icon="bi:info-circle-fill",
                width=22,
            ),
            radius=30,
            size=36,
            variant="outline",
            color="gray",
        ),
        label="More information",
    ),
    id="modal-demo-button",
    n_clicks=0,
)

with open(README_MARKDOWN_FILE, "r", encoding="utf-8") as fp_handle:
    MARKDOWN_TEXT = fp_handle.read()

MODAL_INFO = dmc.Modal(
    title=TAB_TITLE,
    id="modal",
    size="full",
    children=[dcc.Markdown(MARKDOWN_TEXT)],
)

CENTER_GROUP = dmc.Center(
    dcc.Link(
        [
            dmc.MediaQuery(
                create_home_link(BAR_TITLE),
                smallerThan="sm",
                styles={"display": "none"},
            )
        ],
        href="/",
        style={"paddingTop": 5, "textDecoration": "none"},
    ),
)

RIGHT_GROUP = dmc.Group(
    position="right",
    align="center",
    spacing="xl",
    children=[
        GIT,
        MORE_INFO,
        MODAL_INFO,
    ],
)

DASHBOARD_HEADER = dmc.Header(
    height=70,
    p="md",
    children=[
        dmc.Container(
            fluid=True,
            children=dmc.Group(
                position="apart",
                align="flex-start",
                children=[
                    CENTER_GROUP,
                    RIGHT_GROUP,
                ],
            ),
        )
    ],
    mb=25,
)


@app.callback(
    Output("modal", "opened"),
    Input("modal-demo-button", "n_clicks"),
    State("modal", "opened"),
    prevent_initial_call=True,
)
def modal_demo(_, opened):
    """Open README"""
    return not opened
