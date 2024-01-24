from Social.PoliticsWW.dash import Dash, html

app = Dash(__name__)

app.layout = html.Div([
    html.Div(children='Political Regimes')
])

if __name__ == '__main__':
    app.run(debug=True)
