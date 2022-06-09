 function pdbe_get() {
 
    Param([string]$pdb_code)
    
    $pdb_code = $pdb_code.ToLower()
 
    $uri = [System.Uri]"http://www.ebi.ac.uk/pdbe/entry-files/download/pdb${pdb_code}.ent"

    $pdb_file = Join-Path -Path $(Convert-Path '.') -ChildPath $uri.Segments[-1]  # Absolute path for output file
 
    $wc = New-Object System.Net.WebClient
 
    $wc.DownloadFile($uri, $pdb_file)
}