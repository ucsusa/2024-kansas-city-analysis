# PDF to Table

This document describes how I extracted a table of community area boundaries from
[this PDF](https://www.wycokck.org/files/assets/public/v/1/neighborhood-resource-center/documents/2011-neighborhood-directory.pdf) from the Livable Neighborhoods Task Force in Wyandotte County, KS.


## Attempt 1

Following the simple tutorial on [Geeksforgeeks](https://www.geeksforgeeks.org/how-to-extract-pdf-tables-in-python/#) I tried the `tabula-py` and `tabulate` libraries.

```bash
python -m pip install tabula-py tabulate
```

```python
df = read_pdf("path_to_doc.pdf", pages="9-15")
```

`tabula-py` parsed this as multiple tables (one for each page), the results improved when I specified

```python
df = read_pf("path_to_doc.pdf", pages="9-15", multiple_tables=False, output_format="dataframe")[0]
display(df)
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Neighborhood Group NORTH</th>
      <th>Unnamed: 1</th>
      <th>SOUTH</th>
      <th>EAST</th>
      <th>Unnamed: 4</th>
      <th>WEST</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Active Citizens Task Force Haskell</td>
      <td>NaN</td>
      <td>Parallel 27th</td>
      <td>NaN</td>
      <td>34th</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>American Heritage Assoc. Independence</td>
      <td>NaN</td>
      <td>Wagon Trail 109th</td>
      <td>NaN</td>
      <td>111th</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Argentine Betterment Corp.</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Argentine Neighborhood Railroad</td>
      <td>NaN</td>
      <td>Steele Rd - Ruby 18th St. - 7th</td>
      <td>NaN</td>
      <td>I635 - 18th</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Development Assoc.</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>228</th>
      <td>of 113th</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>229</th>
      <td>Wyandotte Gardens Shawnee Drive</td>
      <td>NaN</td>
      <td>Woodend 27th</td>
      <td>NaN</td>
      <td>34th</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>230</th>
      <td>Wyandotte Countians Minnesota Ave.</td>
      <td>NaN</td>
      <td>I-70 10th</td>
      <td>NaN</td>
      <td>40th</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>231</th>
      <td>Against Crime</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>232</th>
      <td>Wyandotte Village Parallel</td>
      <td>NaN</td>
      <td>State N. 42nd</td>
      <td>NaN</td>
      <td>N. 47th</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
<p>233 rows × 6 columns</p>
</div>

Clearly, `tabula-py` is having an issue distinguishing the different columns, even though the rows appear to make sense.

## Attempt 2: `camelot`

After some googling for alternatives, I found the Python package [`camelot-py`](https://camelot-py.readthedocs.io/en/master/). After checking the "comparisons" page that showed performance benchmarks from other pdf reading tools, I tried the following.

```python
import camelot
```
```python
ImportError: DLL load failed while importing cv2. The specified module could not be found.
```

Googling the error message and finding [this post](https://stackoverflow.com/questions/43184887/dll-load-failed-error-when-importing-cv2) on Stackoverflow, I installed

```bash
python -m pip install opencv-contrib-python
```

which resolved the import issue. Now

```py
import camelot
tables = camelot.read_pdf("path_to_doc.pdf")
display(tables[0].df)
```
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
      <th>1</th>
      <th>2</th>
      <th>3</th>
      <th>4</th>
      <th>5</th>
      <th>6</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>N e i g h b o r h o o d s</td>
      <td></td>
      <td></td>
      <td>.\nt o … .\nW e l c o m e  \nSummer 2011 \nt e...</td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>1</th>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>2</th>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>3</th>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>4</th>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>5</th>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>6</th>
      <td></td>
      <td></td>
      <td>Funded by: Community Development Block Grant F...</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
  </tbody>
</table>
</div>

This is extremely wrong. But I think I can correct the issue with fiddling.

## Attempt 3: `camelot` continued

### Specifying page numbers

Next, I specified page numbers

```py
import camelot
tables = camelot.read_pdf("path_to_doc.pdf", pages="9-15")
display(tables)

# <TableList n=10>

display(tables[0].df.head())
```
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
      <th>1</th>
      <th>2</th>
      <th>3</th>
      <th>4</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Neighborhood Group</td>
      <td>NORTH</td>
      <td>SOUTH</td>
      <td>EAST</td>
      <td>WEST</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Active Citizens Task Force  Haskell</td>
      <td></td>
      <td>Parallel</td>
      <td>27th</td>
      <td>34th</td>
    </tr>
    <tr>
      <th>2</th>
      <td>American Heritage Assoc.</td>
      <td>Independence</td>
      <td>Wagon Trail</td>
      <td>109th</td>
      <td>111th</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Argentine Betterment Corp.</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>4</th>
      <td>Argentine Neighborhood \nDevelopment Assoc.</td>
      <td>Railroad</td>
      <td>Steele Rd - Ruby</td>
      <td>18th St. - 7th</td>
      <td>I635 - 18th</td>
    </tr>
  </tbody>
</table>
</div>

That's a huge improvement over the first two attempts!
* `camelot` understood that there were only 5 columns and,
* mostly got the values for each column correctly!

Except, it returned 10 different tables, not all of which were relevant.
Additionally, the four "header" columns we saw in the first table was not
preserved. More fiddling is needed!


### Reading the advanced usage documentation

Since the requirements of this parsing seem to be somewhat more involved, I 
read through the [advanced usage documentation](https://camelot-py.readthedocs.io/en/master/user/advanced.html#specify-column-separators) and tried a few things to learn more about
how the parser was interpreting the document. There were some very interesting plots!

Some notes:

* The "flavors" `lattice` and `stream` produce slightly different results. The former uses an image of the pdf to parse, while the latter tries a line-by-line approach (i.e., a stream).


Ultimately, I landed on this

```py
flavor='stream'

x1 = 25
y1 = 670
x2 = 580
y2 = 75
cols = [160, 275, 375, 460]
p1 = 9
p2 = 15
tables = camelot.read_pdf(str(pdf_path), 
                          pages=f"{p1}-{p2}", 
                          flavor=flavor,
                          table_areas=[f'{x1},{y1},{x2},{y2}'],
                          columns=[','.join([str(col) for col in cols])],
                          row_tol=15,
                          column_tol=0,
                          strip_text='.\n',
                          )
```