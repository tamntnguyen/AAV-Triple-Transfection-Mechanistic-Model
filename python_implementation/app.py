#!/usr/bin/env python

"""
 Creates a Data Science Dashboard for
 the AAV ODE model.

 + Uses Plotly Dash (https://dash.plotly.com)
 + 2 tabs:  Data tab and Figures tab

Usage:
 python app.py

"""
from pandas.core.frame import DataFrame

from dash import html, dash_table, dcc, callback_context
from dash.dependencies import Input, Output

import dash_mantine_components as dmc
from dash_iconify import DashIconify

from header_bar import app, DASHBOARD_HEADER

from aav_ode import Experiments

aav_model = Experiments()

MAX_HOURS = 120
default_media_change = aav_model.get_media_exchange_time()
default_total_time = aav_model.get_total_time()
initial_concentrations = aav_model.get_initial_concentrations()


################################
# PUBLIC FUNCTIONS
################################


def show_data(input_value):
    """Show the data tables

    Creates a webpage with the data table(s)

    Args:
        input_value(int): The value of the LED slider (outgrowth threshold)

    Returns:
        A Dash table
    """
    aav_model.set_total_time(input_value[1])
    aav_model.set_media_exchange_time(input_value[0])
    time_pt, prediction = aav_model.run_simulation()

    df_data = DataFrame(prediction.T, columns=aav_model.get_output_labels())
    df_data.insert(0, "Hours post transfection", time_pt)
    dt_col_param = []

    for col in df_data.columns:
        dt_col_param.append(
            {
                "name": str(col),
                "id": str(col),
                "type": "numeric",
                "format": dash_table.Format.Format(precision=2),
            }
        )

    table = dash_table.DataTable(
        columns=dt_col_param,
        data=df_data.to_dict("records"),
        sort_action="native",
        filter_action="native",
        page_action="native",
        page_current=0,
        page_size=20,
        style_table={"overflowX": "auto", "overflowY": "auto"},
        export_format="xlsx",
    )
    return dmc.Paper(children=[table])


def show_figures(input_value):
    """Show the figures

    Creates a webpage with the graph(s)


    Returns:
        A list of Dash sections
    """
    aav_model.set_total_time(input_value[1])
    aav_model.set_media_exchange_time(input_value[0])

    time_pt, prediction = aav_model.run_simulation()

    webpage = [
        dmc.Paper(
            children=[
                dmc.Paper(
                    children=dcc.Graph(
                        id="replication",
                        figure=aav_model.plot_replication(time_pt, prediction),
                    )
                ),
                dmc.Paper(
                    children=dcc.Graph(
                        id="production",
                        figure=aav_model.plot_production(time_pt, prediction),
                    ),
                ),
            ]
        ),
        dmc.Paper(
            children=[
                html.Div("Components"),
                dcc.Graph(
                    id="plot_outputs",
                    figure=aav_model.plot_outputs(time_pt, prediction),
                ),
            ],
        ),
    ]

    return webpage


################################
# Dash Plotly Webpage Layout
################################

range_slider = dcc.RangeSlider(
    id="time-window-slider",
    min=0,
    max=MAX_HOURS,
    step=0.5,
    allowCross=False,
    marks={
        int(default_media_change): {
            "label": "Default media exchange",
            "style": {"color": "green"},
        },
        int(default_total_time): {
            "label": "Default total time",
            "style": {"color": "green"},
        },
        MAX_HOURS: {
            "label": "Maximum hours",
            "style": {"color": "red"},
        },
    },
    tooltip={"placement": "top", "always_visible": True},
    value=[
        aav_model.get_media_exchange_time(),
        aav_model.get_total_time(),
    ],
)

PAPER_URL = (
    "https://www.cell.com/"
    "molecular-therapy-family/methods/fulltext/"
    "S2329-0501(21)00072-3"
)

"""
Add model settings (k constants) to settings page
Allow user to modify the k settings.
"""


def generate_setting_field(input_key, input_value):
    """
    Create input widgets for the settings
    """

    val = aav_model.get_k(input_key)

    return dmc.SimpleGrid(
        cols=2,
        children=[
            dmc.Text(input_value),
            dcc.Input(value=val, min=1e-10, max=1e10, type="number", id=input_value),
        ],
    )


settings_page = dmc.Container(
    dmc.Group(
        children=[
            generate_setting_field(key, value)
            for key, value in aav_model.get_input_labels().items()
        ],
        position="left",
    )
)

offcanvas = dmc.Drawer(
    children=[settings_page],
    id="offcanvas",
    title="AAV Model Settings",
    size="lg",
    padding="md",
    shadow="md",
)

settings_button = dmc.Button(
    "Settings",
    id="open-offcanvas",
    variant="gradient",
    gradient={"from": "teal", "to": "blue", "deg": 60},
    leftIcon=[DashIconify(icon="fluent:settings-32-regular")],
)

initial_conditions = dmc.Group(
    [
        dmc.NumberInput(
            label="Initial Concentration (plasmid/cell)",
            description="From 0 to infinity",
            value=initial_concentrations[0],
            min=0,
            max=1e6,
            id="initial0",
        ),
        dmc.NumberInput(
            label="Initial Concentration (plasmid/cell)",
            description="From 0 to infinity",
            value=initial_concentrations[1],
            min=0,
            max=1e6,
            id="initial1",
        ),
        dmc.NumberInput(
            label="Initial Concentration (plasmid/cell)",
            description="From 0 to infinity",
            value=initial_concentrations[2],
            min=0,
            max=1e6,
            id="initial2",
        ),
    ]
)


figures_and_data = dmc.Tabs(
    [
        dmc.Tab(label="Figures", id="figures"),
        dmc.Tab(label="Data", id="data"),
    ],
    id="tabs",
    active=0,
)

# Create an input list of everything on the settings page
items = []
for key, value in aav_model.get_input_labels().items():
    items.append(Input(component_id=value, component_property="value"))

items.append(Input("tabs", "active"))
items.append(Input("time-window-slider", "value"))
items.append(Input("initial0", "value"))
items.append(Input("initial1", "value"))
items.append(Input("initial2", "value"))


@app.callback(
    Output("tab-content", "children"),
    [items],
)
def update_settings(*args):
    """Update the model settings

    Triggered whenever a setting is changed by the user.
    """

    n_settings = len(aav_model.get_input_labels().items())
    active_tab = args[n_settings]
    slider_value = args[n_settings + 1]
    initial0 = args[n_settings + 2]
    initial1 = args[n_settings + 3]
    initial2 = args[n_settings + 4]

    aav_model.set_initial_concentrations([initial0, initial1, initial2])
    # Figure out which setting was changed
    ctx = callback_context

    # Get the triggered setting name
    triggered_name = ctx.triggered[0]["prop_id"].split(".")[0]
    new_value = ctx.triggered[0]["value"]

    if new_value is not None:

        idx = 0
        for _, name in aav_model.get_input_labels().items():

            if name == triggered_name:
                # Update setting with new value
                aav_model.set_k(idx, float(new_value))

            idx += 1

    if active_tab == 1:
        return show_data(slider_value)

    if active_tab == 0:
        return show_figures(slider_value)

    raise Exception(f"Tab {active_tab} is undefined here.")


@app.callback(
    Output("offcanvas", "opened"),
    Input("open-offcanvas", "n_clicks"),
    prevent_initial_call=True,
)
def open_drawer(_):
    """Open settings page"""
    return True


layout = dmc.Container(
    [
        DASHBOARD_HEADER,
        range_slider,
        html.Br(),
        dmc.Grid(
            [
                dmc.Col(initial_conditions, span=10),
                dmc.Col(settings_button, span=1),
            ]
        ),
        offcanvas,
        html.Br(),
        figures_and_data,
        html.Div(id="tab-content", className="p-4"),
    ],
    fluid=True,
)

app.layout = dmc.Container(children=layout)

server = app.server  # expose server variable for Procfile

if __name__ == "__main__":
    app.run_server(debug=True)
