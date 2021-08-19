# Analysis of SARS-CoV-2 deep sequencing from Shen et al (2020)

## Background

There are discrepancies in dates in the sequence collection for the eight samples: they were originally described as being from December 18 to 29, 2019.
But then at some post peer-review stage after posting of the early access version, this was changed to January 2020.
Specifically:

### Final published version says January 2020
The final version of [Shen et al (2020)](https://academic.oup.com/cid/article/71/15/713/5780800) was published in _Clinical Infectious Diseases_ on May-5-2020.
That final version (archived on the Wayback Machine [here](https://web.archive.org/web/20210628200843/https://academic.oup.com/cid/article/71/15/713/5780800)) describes the sequences coming from _"Eight COVID-19 pneumonia samples were collected from hospitals in Wuhan on January 2020"_.

### PubMed Central version says December 18 to 29, 2019

However, the [PubMed Central version of Shen et al (2020)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7108196/pdf/ciaa203.pdf) describes the samples as _"Eight COVID-19 pneumonia samples were collected from hospitals in Wuhan from December 18 to 29, 2019"_. This PubMed Central version appears to have been built from the _Clinical Infectious Diseases_ original online-early access version that posted on March-9-2020.
I archive the Pubmed Central version on the Wayback Machine [here](https://web.archive.org/web/20210817201634/https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7108196/pdf/ciaa203.pdf), and also [downloaded a copy of the PDF here](data/ciaa203.pdf).

### Early access version from journal says December 18 to 29, 2019

The overall article history provided by _Clinical Infectious Diseases_ indicates the paper was received Feb-18-2020, accepted on Feb-25-2020, first published on March-9-2020, and corrected and typeset on May-5-2020.
I was able to find the early access version for _Clinical Infectious Diseases_ archived on the Wayback Machine [here](https://web.archive.org/web/20200308003758/https://watermark.silverchair.com/ciaa203.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAmswggJnBgkqhkiG9w0BBwagggJYMIICVAIBADCCAk0GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMUvB1YKwP2bcB5dTjAgEQgIICHhCVB_L3-QgfUyzQ5G-v68bJ1kHhDf8lLShRGqC9bK7DlebkU5yUv6DPXGfkYCHNhvUVkawokppu_0H1qyLMACdOSYBspS7NvV5gembGXwlwt9_Aci4waeoxJX96d777yav5KH7giT_Tx8xzQR26RlMHyFDq4oUhCej7jiTdl0H3JFPULiDU1FGeW9SA3EiSy4iKiLv4n2lyTBkTYF2LpQ59SKThu1vic0sPDTpTihihfPedjRnk3XoIHIWWRzO_SmC7gBYZs6i5kY-XfDUBRauvezNPPOcrwC7bcUCjJ4OEKa8ym9ssdhVlKHgB-UPbQLRaTZ285BTc5f1wWAjRRJX7Hkb7ik4A5qWkKnFePKsu_rRuJ_8cKRv9rfYJODY2VNeHGf3QyHrIC_xCSMjPRlTGpTD5aTdw2vIQuPzf9hLT44YBGKD0beExJUhqCMeq2i5yUn582z8MMW4LYKq-0qS49z0mjAGVjvpftdPVY-6MhhjL7KWz48RB69vFaMTqybNKmqKOg0NKxXLi4iq8C3orPY_sBKagy1GKGuvfz39pUkh6MyAwM1XPzKmzevfDKd8xHNyYGdRTK8F3ThZI4PTdvkMXSq0Tp8v3d3x6mpsUzQ0bd2ypmQ7KyUauEnH_AxfZ8nYIYq9KebygrgFHVqkHzUb2l53vE0rbYqQdy1twbFJCkIfGtRpIkQfl5izggBvLKWOw5yA7sLDrJJba), and it also refers to the samples as from _"December 18 to 29, 2019"_.

The supplementary data links for the early access version is archived [here](https://web.archive.org/web/20200308003745/https://academic.oup.com/cid/advance-article/doi/10.1093/cid/ciaa203/5780800?searchresult=1#supplementary-data).
However, the actual links to the supplementary material are broken, so I cannot access the actual supplementary data for the early access version.

## Analysis pipeline
The analysis pipeline in [Snakefile](Snakefile) gets the accessions from both the SRA uploaded BAMs and SRA files, and alignms them to Wuhan-Hu-1 and calls variants.
Before doing this, it splits into separate files / Illumina runs and analyzes those independently.

The final variant calls are in [results/aggregated_variants/all_samples.csv](results/aggregated_variants/all_samples.csv).

Then the Jupyter notebook [analysis.ipynb](analysis.ipynb) is run manually outside the `snakemake` pipeline to look at mutations.
Although some of the samples have a lot of mutations, they are confusing and don't seem to be ones to more ancestral variants (nothing at 8782, 18060, 28144, or 29095).
