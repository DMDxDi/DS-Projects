import pandas as pd
from bs4 import BeautifulSoup
import requests

html_text = requests.get('https://en.wikipedia.org/wiki/Polity_data_series').text

soup = BeautifulSoup(html_text, 'lxml')

countries = []
regimes = []

# Find the table with the data
table = soup.find('table', {'class': 'wikitable'})

# Find all rows in the table
rows = table.find_all('tr')

for row in rows[1:]:  
    columns = row.find_all('td')
    
    if len(columns) >= 6: 
        country = columns[0].get_text(strip=True)
        regime = columns[5].get_text(strip=True)
        
        countries.append(country)
        regimes.append(regime)

# Save data
data = {'Country Name': countries, 'Regime': regimes}
df = pd.DataFrame(data)
df.to_csv('political_regimes.csv', index=False)

