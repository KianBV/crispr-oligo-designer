# Microsynth template folder

Place the official Microsynth DNA order-form template here:

```
templates/MicrosynthUploadFormDNA__2_.xlsx
```

The application fills a **copy** of this template — it never modifies the
original. The path is set in [`config/config.yml`](../config/config.yml) under
`microsynth.template_path`. You can also upload a template at runtime from the
app's **Settings** panel, which overrides the configured path for that session.

## Expected layout (sheet "DNA Order")

The app writes into these columns, starting at row 2 (row 1 is the header):

| Col | Field                | Written by the app                     |
|-----|----------------------|----------------------------------------|
| A   | Name                 | `<mutation>_pEcgRNA_fw_gibson/_rv_gibson/_HA` |
| B   | Sequence             | oligo sequence (PTO `*` on homology arm)|
| C   | Length               | left as the template's formula          |
| D   | Purification         | default `DES`                           |
| E   | Scale                | default `0.04`                          |
| F–K | Modifications        | left empty                              |
| L   | Physical Condition   | default `Dried`                         |
| M   | Datasheet            | default `Standard`                      |
| N   | Aliquots             | default `No`                            |

The column defaults (D, E, L, M, N) are configurable in `config/config.yml`.
