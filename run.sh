#Rscript lib/create_plotly.R ../Apr26-2022/Results/mapped_gc_content_threegenomesRef.csv Apr26 reference
#Rscript lib/create_plotly.R ../Aug3-2022/Results/mapped_gc_content_threegenomesRef.csv Aug3 reference
#Rscript lib/create_plotly.R ../Jul26-2022/Results/mapped_gc_content_threegenomesRef.csv Jul26 reference
Rscript lib/create_plotly.R ../Apr26-2022/Results/mapped_gc_content_threegenomesRef.csv Apr26 Ref
Rscript lib/create_plotly.R ../Aug3-2022/Results/mapped_gc_content_threegenomesRef.csv Aug3 Ref
Rscript lib/create_plotly.R ../Jul26-2022/Results/mapped_gc_content_threegenomesRef.csv Jul26 Ref

Rscript lib/create_plotly.R ../Apr26-2022/Results/mapped_gc_content_threegenomesstrain.csv Apr26 strain
Rscript lib/create_plotly.R ../Aug3-2022/Results/mapped_gc_content_threegenomesstrain.csv Aug3 strain
Rscript lib/create_plotly.R ../Jul26-2022/Results/mapped_gc_content_threegenomesstrain.csv Jul26 strain

tar -czvf plotly.tar.gz ../Aug3-2022/Results/*.html ../Jul26-2022/Results/*.html ../Apr26-2022/Results/*.html
