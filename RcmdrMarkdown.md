<!-- R Commander Markdown Template -->

Replace with Main Title
=======================

### Your Name

### 2020-02-14







```r
> .NewData <- data.frame(age=60, lcavol=1.5, row.names="1")
> .NewData  # Newdata
```

```
  age lcavol
1  60    1.5
```

```r
> predict(RegModel.1, newdata=.NewData, interval="confidence", level=.95, 
+   se.fit=FALSE)
```

```
       fit      lwr     upr
1 2.583387 2.398925 2.76785
```


