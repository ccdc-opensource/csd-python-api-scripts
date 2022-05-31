# Run a docking via the CLI for comparison with those run via the API.
# gold.conf here should be the same as e.g. ..\output_foreground\api_gold.conf except for file paths.

Remove-Item -Recurse -Force -ErrorAction SilentlyContinue .\output

& 'C:\Program Files\CCDC\Discovery_2021\GOLD\gold\d_win32\bin\gold_win32.exe' .\gold.conf

# To start Hermes and view results...

# & 'C:\Program Files\CCDC\Discovery_2021\Hermes\hermes.exe' .\gold.conf