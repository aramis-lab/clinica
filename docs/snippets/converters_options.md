## Optional parameters common to all converters

- `-- subjects_list` / `- sl`  : path to a text file containing a list of specific subjects to extract. The expected format is one subject per line :

=== "Example with the ADNI dataset"
```
  001_S_0001
  002_S_0002
```

- `-- n_procs` / `- np` :  Number of cores used to run in parallel.  (default: (Number of available CPU minus one))
