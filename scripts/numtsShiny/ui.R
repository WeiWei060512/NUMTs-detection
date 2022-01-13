library(shiny)
library(shinydashboard)
library(DT)
library(shinyjs)
library(sodium)
library(data.table)
library(openxlsx)
library(shinyWidgets)

#https://www.listendata.com/2019/06/how-to-add-login-page-in-shiny-r.html

project_name = 'Nuclear-embedded mitochondrial DNA sequences'
#source('module_login.R')
germline_df = as.data.table(read.xlsx('germlineNUMTs.xlsx'))
cancer_df = as.data.table(read.xlsx('cancerNUMTs.xlsx'))

### define lists ###

projects_list = c('Germline', 'Cancer')
population_list = c('African','American','EastAsian','SouthAsian','European','Others')
frequency_list = c('common','rare','ultra-rare')

chromosome_list = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
                    'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
                    'chr18','chr19','chr20','chr21','chr22','chrX')
mtGene_list =  c('Dloop','rRNAs','tRNAs','ATP6','ATP8','CO1','CO2','CO3','CYB','ND1','ND2','ND3','ND4/ND4L','ND5','ND6')
NUMTtype_list = c('UN','validated','putative')
longRead_list = c('Y','Na')
cancerType_list = unique(cancer_df$Tumour_Type)
cancerGene_list = c('Y','N')
complexNUMTs_list = c('UN','complexNUMTs')

b64_germline <- base64enc::dataURI(file="www/germline.sample.png", mime="image/png")
b64_cancer <- base64enc::dataURI(file="www/cancer.sample.png", mime="image/png")

header <- dashboardHeader( title = "NUMTs", uiOutput("logoutbtn"))
sidebar <- dashboardSidebar(uiOutput("sidebarpanel")) 
body <- dashboardBody(shinyjs::useShinyjs(), uiOutput("body"))
dashboardPage(header, sidebar, body, skin = "blue")

