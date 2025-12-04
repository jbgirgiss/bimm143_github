# Class 6: R functions
Joseph Girgiss (PID: A17388247)

All functions in R have at least 3 things:

- A **name**, we pick this and use it to call the function.
- Input **arguments**, there can be multiple comma separated inputs to
  the function.
- The **body**, lines of R code that do the work of the function.

Our first function:

``` r
add <- function(x, y=1) {
  x + y
}
```

Let’s test out function:

``` r
add(c(1, 2, 3), y=10)
```

    [1] 11 12 13

``` r
add(10)
```

    [1] 11

``` r
add(10, 100)
```

    [1] 110

## A second function

Let’s try something more interesting. Make a sequence generation tool.

The `sample()` function can be useful here to work with nucleotides A,
C, T, G and return 3 of them.

``` r
sample(c("A", "T", "C", "G"), size=3)
```

    [1] "T" "C" "G"

``` r
sample(1:10, size=3)
```

    [1]  7 10  3

``` r
n <- c("A", "T", "C", "G")
       sample (n, size=15, replace = TRUE) 
```

     [1] "T" "A" "A" "T" "A" "C" "C" "C" "A" "A" "T" "G" "C" "A" "C"

Turn this snippet into a function that returns a user specified length
DNA sequence. Let’s call it `generate_dna()`…

``` r
generate_dna <- function(length) {
  n <- c("A", "T", "C", "G")
  v <- sample(n, size=length, replace = TRUE)
  cat("Well done you!\n")
  return(v)
  }
```

``` r
generate_dna(30)
```

    Well done you!

     [1] "T" "C" "A" "A" "G" "A" "C" "A" "A" "G" "G" "A" "T" "A" "C" "A" "T" "C" "A"
    [20] "A" "T" "C" "G" "T" "T" "A" "C" "C" "T" "C"

``` r
s <- generate_dna(15)
```

    Well done you!

``` r
s
```

     [1] "T" "G" "C" "T" "G" "T" "A" "T" "G" "A" "G" "G" "C" "G" "A"

I want the option to return a single element character vector with my
sequence all together like this: “GGAGTAC”

``` r
generate_dna <- function(length, fasta=FALSE) {
  n <- c("A", "T", "C", "G")
  v <- sample(n, size=length, replace = TRUE)
  
  # Make a single element vector
  s <- paste(v, collapse = "")
  
  cat("Well done you!\n")
  
  if (fasta) {
    return(s)   # return s if fasta=TRUE
  } else {
    return(v)   # return v if fasta=FALSE
  }
}
```

``` r
generate_dna(10, fasta=TRUE)
```

    Well done you!

    [1] "CGCTTAGTGT"

## A more advanced example

Make a third function that generates protein sequence of a user
specified length and format.

``` r
generate_protein <- function(length, collapse = TRUE) {
  # Define the 20 standard amino acids
  amino_acids <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", 
                   "I", "L", "K", "M", "F", "P", "S", "T", "W", 
                   "Y", "V")
  
  # Generate a random sequence
  seq <- sample(amino_acids, size = length, replace = TRUE)
  
  # Return as string or vector depending on user input
  if (collapse) {
    return(paste(seq, collapse = ""))
  } else {
    return(seq)
  }
}
```

``` r
# Generate a 15-amino-acid protein sequence as a string
generate_protein(15)
```

    [1] "QLQLYYLRAPMNVTW"

``` r
# Generate a 10-amino-acid protein sequence as a vector
generate_protein(10, collapse = FALSE)
```

     [1] "V" "E" "L" "P" "W" "Y" "I" "G" "K" "V"

``` r
generate_protein()
```

    [1] "CVYSFHAGQPHIVVWTVQNK"

> Q. Generate a random protein sequences between lengths 5 and 12
> amino-acids.

One approach to do this is by brute force calling our function for each
length 5 to 12.

Another approach is to write a `for()` loop to itterate over the input
valued 5 to 12.

A very useful third R specific approach is to use the `sapply()`
function.

``` r
seq_lengths <- 5:12 
for (i in seq_lengths) {
  cat(">", i, "\n")
  cat(generate_protein(i))
  cat("\n")
}
```

    > 5 
    RVMLN
    > 6 
    GDFQME
    > 7 
    FINDTSP
    > 8 
    WNNLDRTA
    > 9 
    VGIHEIRLQ
    > 10 
    YYQDQCHTEF
    > 11 
    FHWSMHYDHVH
    > 12 
    SLMQKFQHQWCV

``` r
sapply(5:12, generate_protein)
```

    [1] "AWMMP"        "IAWYAC"       "LEICAWN"      "HSRWITYW"     "RHYVTPQWS"   
    [6] "NDKFWLTLIM"   "VVQETLNPPWD"  "AIEMGMKPAVAI"

> **Key-Point**: Writing functions in R is doable but not the easiest
> thing. Start with a working snippet of code and then use LLM tools to
> improve and generalize your function code is a productive approach.
