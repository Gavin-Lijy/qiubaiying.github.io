[toc]

# Introduction

## The challenges to standardization

- Growing number of
  analysis methods and exploding
  dataset sizes
- Analysis tools for scRNA-seq data are written in a variety of
  programming languages

Instead of targeting a standardized analysis pipeline, we outline current best practices
and common tools independent of programming language 



##  Key elements of an experimental scRNA-seq workflow

*Ziegenhain et al (2017);Macosko et al (2015); Svensson et al (2017).*



![image-20200717161034833](https://i.loli.net/2020/07/17/wGbn12FrEQxjNJK.png)

<center>figure 1</center>



### single-cell isolation



#### methods

1. plate-based techniques isolate cells into wells on a plate
2. droplet-based methods rely on capturing each cell in its own microfluidic droplet

#### error

errors can occur that lead to
multiple cells being captured together (doublets or multiplets), nonviable
cells being captured, or no cell being captured at all (empty
droplets/wells).

### library construction

(Library construction is the process in which the intracellular mRNA is captured, reverse-transcribed to cDNA molecules and  amplified)

1. Each well or droplet contains the necessary chemicals to break down the cell membranes and perform library construction

2. the mRNA from each cell can
   be labelled with a well- or droplet-specific cellular barcode

3. experimental protocols also label captured molecules (mRNA) with a unique molecular identifier (UMI).

   (**UMIs allow us to distinguish between amplified copies of the same mRNAmolecule and reads from separate mRNA molecules transcribed from the same gene.**)


### sequencing

#### Without UMI

Sequencing produces **read data**, which undergo quality control, grouping based on their assigned barcodes (demultiplexing) and alignment in read processing pipelines.

#### With UMI

For UMI-based protocols, read data can be further demultiplexed to produce counts of captured mRNA molecules
(**count data**).

# Pre-processing and visualization

- Tools

Cell Ranger (Zheng et al, 2017), indrops (Klein et al, 2015),
SEQC (Azizi et al, 2018), or zUMIs (Parekh et al, 2018)

- Procedure

read quality control (QC)

demultiplexing

（assigning reads to their cellular barcodes and mRNA molecules of origin ）

genome alignment

quantification

## 1. Quality control

![image-20200710150515627](https://i.loli.net/2020/07/10/GolXmDpF4YyqibN.png)

<center>figure2</center>

### Goal

 Ensuring that all cellular barcode data correspond to **viable** cells.

### Procedure

#### Checking  the integrity of cells

There QC covariates（jointly）

1. the number of counts per barcode (count depth)
2.  the number of genes per barcode
3.  the fraction of counts from mitochondrial genes per barcode

##### Tools

*DoubletDecon: preprint: DePasquale et al, 2018; Scrublet: Wolock et al, 2019;Doublet Finder: McGinnis et al, 2018*

#### Checking the level of transcripts

<font color=red>A guideline to setting this threshold is to use the minimum cell cluster size that is of interest and leaving some leeway for dropout effects</font>

#### Correcting the ambient gene expression

**Ambient gene expression** refers to counts that do not originate from a barcoded cell, but from other lysed cells whose mRNA contaminated the cell suspension prior to library construction.

##### Tools

*SoupX (preprint: Young & Behjati, 2018)*

*Angelidis et al, 2019*

### Pitfalls & recommendations

-  Perform QC by finding outlier peaks in the number of genes, the count depth and the fraction of mitochondrial reads. Consider these covariates jointly instead of separately.
-  Be as permissive of QC thresholding as possible, and revisit QC if downstream clustering cannot be interpreted.
-  If the distribution of QC covariates differ between samples, QC thresholds should be determined separately for each sample to account for sample quality differences as in Plasschaert et al (2018).

## 2. Normalization

### Tools

- bulk gene expression

  preprint: Pachter, 2011; Dillies et al, 2013

- single-cell data

  Lun et al, 2016a; Vallejos et al, 2017

### Procedure

####  Cellular count data

##### Count depth scaling (CPM)

  - also referred to as “counts per million” or CPM normalization.
  - This protocol comes from bulk expression analysis and normalizes count data using a so-called **size factor** proportional to the count depth per cell.
  - assuming that all cells in the dataset initially contained an equal number of mRNA molecules and count depth differences arise only due to sampling.

##### Extension of CPM

- high-count filtering CPM

  *Weinreb et al (2018)*

- *Scran’s pooling-based size factor estimation method (Lun et al, 2016a)*

  size factors are estimated based on a **linear regression** over genes, after cells are pooled to avoid technical dropout effects.

##### Non-linear normalization methods

- particularly relevant for **plate-based** scRNA-seq data, which tend to have batch effects between plates.

- Such an approach can combine technical and biological data correction (e.g. batch correction or correction for cell cycle effects) with count depth normalization.

###### Tools

- Mayer et al (2018) 

  it a negative binomial model to count data, using technical covariates such as the read depth and the number of counts per gene to fit the model parameters.The residuals of the model fit serve as a normalized quantification of gene expression.

- *Cole et al, 2019*

  outperform global scaling methods especially in situations with strong batch effects 

##### Appropriate normalization method for a specific dataset

- *Cole et al (2019)*

  their **scone** tool should be used to select the appropriate normalization method for a specific dataset

- ScRNA-seq techniques can be divided into **full-length** and **3' enrichment** methods (Svensson et al, 2017; Ziegenhain et al, 2017).

###### Full-length

- benefit from normalization methods that take into account gene length (e.g. Patel et al, 2014; Kowalczyk et al, 2015; Soneson & Robinson, 2018)

- TPM normalization (Li et al, 2009)

###### 3' enrichment

scran for normalization of non-full-length datasets

#### Gene counts

- Gene normalization constitutes scaling gene counts to have zero mean and unit variance (z scores)

- There is currently no consensus on whether or not to perform normalization over genes.

1. pro:

Seurat tutorials (Butler et al, 2018) generally
apply gene scaling

2. con 

the authors of the Slingshot method opt against
scaling over genes in their tutorial (Street et al, 2018).

- The preference between the two choices revolves around whether all genes should be weighted equally for downstream analysis, or whether the magnitude of expression of a gene is an informative proxy for the importance of the gene.
- opt to refrain from scaling over genes in this tutorial.

#### log(x+1) transformation

##### advantage

- effect

1. distances between log-transformed expression values represent log fold changes, which are the canonical way to measure changes in expression.
2. log transformation mitigates (but does not remove) the mean–variance relationship in single-cell data (Brennecke et al, 2013)
3. log transformation reduces the skewness of the data to approximate the assumption of many downstream analysis tools that the data are **normally distributed.**

- This usefulness is highlighted by down stream applications for **differential expression testing** (Finak et al, 2015; Ritchie et al, 2015) or **batch correction** (Johnson et al, 2006; Buttner et al, 2019) that use

##### disadvantage

introduce spurious differential expression effects into the
data (preprint: Lun, 2018). 

This effect is particularly pronouncedwhen normalization size factor distributions differ strongly between tested groups.

### Pitfalls & recommendations:



- We recommend scran for normalization of non-full-length datasets. An alternative is to evaluate normalization approaches via scone especially for plate-based datasets. 
- Full-length scRNA-seq protocols can be corrected for gene length using bulk methods.
-  Normalized data should be log(x+1)-transformed for use with downstream analysis methods that assume data are normally distributed.

## 3. Data correction and integration

### Goal

- Data correction targets further technica and biological covariates such as **batch, dropout, or cell cycle effects.**

- the decision of which covariates to consider will depend on the intended downstream analysis.

### Methods

#### Regressing out biological effects

##### Remove the effects of the cell cycle on the transcriptome

- tools

implemented in the Scanpy and Seurat platforms (Butler
et al, 2018; Wolf et al, 2018) 

 in specialized packages with more complex mixture models such as scLVM (Buettner et al, 2015) or fscLVM (Buettner et al, 2017).

- also be used to regress out other known biological
  effects such as mitochondrial gene expression, which is interpreted as an indication of cell stress.

##### Aspects considered 

- cell cycle signals can also be informative of the biology

- correcting for one process may unintentionally mask the signal of another.
- correcting for cell size via normalization, or dedicated tools such as cgCorrect (Blasi et al, 2017), also partially corrects for cell cycle effects in scRNA-seq data.

#### Regressing out technical effects

- Most prominent technical covariates are **count depth and batch**
- Regressing out count depth effects can improve the performance of trajectory inference algorithms, which rely on finding transitions between cells

- An alternative to regression-based strategies for removing count effects is to use a more rigorous normalization procedure such as downsampling or non-linear normalization methods

#### Batch effects and data intergration

- Three levels of batch effects

1. between groups of cells in an experiment,

**batch correction**  linear methods

2.between experiments performed in the same laboratory or between datasets from different laboratories.

**data integration. ** non-linear approaches

##### Batch effects

![image-20200717161243984](https://i.loli.net/2020/07/17/XyL8jUDNAsvCGrZ.png)

<center>figure 3 </center>



- *tools: ComBat (Johnson et al, 2006)*

ComBat uses all cells in a batch to fit batch parameters. This approach will confound the batch effect with biological differences between cell types or states that are not shared among datasets

- experimental methods

best method of batch correction is pre-empting the effect and avoiding it altogether by clever experimental design (Hicks et al,2017) Batch effects can be avoided by pooling cells across experimental conditions and samples. Using strategies such as cell tagging (preprint: Gehring et al, 2018), or via genetic variation (Kang et al, 2018), it is possible to demultiplex cells that were pooled in the experiment.

##### Data integration

- tools

*Canonical Correlation Analysis (CCA; Butler et al, 2018), Mutual Nearest Neighbours (MNN; Haghverdi et al, 2018), Scanorama (preprint: Hie et al, 2018), RISC (preprint: Liu et al, 2018), scGen (preprint: Lotfollahi et al, 2018), LIGER (preprint: Welch et al, 2018), BBKNN (preprint: Park et al, 2018), and Harmony (preprint: Korsunsky et al, 2018)*

#### Expression recovery

Also denoising or imputation

- Inferring dropout events, replacing these zeros with appropriate expression values, and reducing the noise in the dataset

- tools

  *(MAGIC: van Dijk et al, 2018; DCA: Eraslan et al, 2019; scVI: Lopez et al, 2018; SAVER: Huang et al, 2018; scImpute: Li & Li, 2018)*

- this step can be integrated with normalization, batch correction and other downstream analysis as implemented in the scVI tool (Lopez et al, 2018).
- There is currently no consensus on how denoised data should be used in the light of these considerations

### Pitfalls & recommendations

- Regress out biological covariates only for trajectory inference and if other biological processes of interest are not masked by the regressed out biological covariate.
- Regress out technical and biological covariates jointly rather than serially.
- We recommend performing batch correction via ComBat when cell type and state compositions between batches are consistent

- Data integration and batch correction should be performed by different methods. Data integration tools may over-correct simple batch effects.
- Users should be cautious of signals found only after expression recovery. Exploratory analysis may be best performed without this step.

## 4. Feature selection, dimensionality reduction and visualization

### Feature selection

#### Goal

keep only genes that are “informative” of the variability in the
data (highly variable genes (HVGs))

#### Tools

highly variable genes (HVGs)Brennecke et al, 2013

There exist different flavours of this algorithm that expect either count data (Seurat) or log-transformed data (Cell Ranger)

reviewed in Yip et al (2018).

### Dimentionally reduction

![image-20200716222136379](https://i.loli.net/2020/07/16/uy2wDltoIO6UVzK.png)

<center> figure 4</center>

#### Goal

Visualization and summarization

#### Tools

##### PCA

*principal component analysis(PCA; Pearson, 1901)* 

we can orrelate quantities of interest with principal components to assess their importance.

##### diffusion maps

*(Coifman et al, 2005)*

As diffusion components emphasize transitions in the data, they are principally used when continuous processes such as differentiation are of interest. Typically, each diffusiocomponent (i.e. diffusion map dimension) highlights the heterogeneity of a different cell population.

### Visualization

#### Tools

##### t-SNE

- t-distributed stochastic neighbour embedding (t-SNE; van derMaaten & Hinton, 2008).

- t-SNE dimensions focus on capturing local similarity at the expense of global structure
- disadvantage
  1. these visualizations may exaggerate differences between cell populations and overlook potential connections between these populations
  2. A further difficulty is the choice of its perplexity parameter, as t-SNE graphs may show strongly different numbers of clusters depending on its value (Wattenberg et al, 2016)

##### UMAP 

- Common alternatives to t-SNE are the Uniform Approximation and Projection method (UMAP; preprint: McInnes & Healy, 2018) or graph-based tools such as
  SPRING (Weinreb et al, 2018).
- force-directed layout algorithm ForceAtlas2 arguably represent the best approximation of the underlying topology (Wolf et al, 2019,
- advantage
  1. What sets UMAP apart in this comparison is its speed and ability to scale to large numbers of cells (Becht et al, 2018).
  2. UMAP can also summarize data in more than two dimensions

##### PAGA

- partition-based graph abstraction (PAGA; Wolf et al, 2019

- adequately approximate the topology of the data while coarse-graining the visualization using clusters. In combination with any of the above visualization methods, PAGA produces coarse-grained visualizations, which can simplify the
  interpretation of single-cell data especially with large numbers of cells.

### Pitfalls & recommendations

- We recommend selecting between 1,000 and 5,000 highly variable genes depending on dataset complexity.
- Feature selection methods that use gene expression means and variances cannot be used when gene expression values have been normalized to zero mean and unit variance, or when residuals from model fitting are used as normalized expression values.
- <font color=red>Dimensionality reduction methods should be considered separately for summarization and visualization.</font>
- We recommend UMAP for exploratory visualization; PCA for general purpose summarization; and diffusion maps as an alternative to PCA for trajectory inference summarization.
- PAGA with UMAP is a suitable alternative to visualize particularly complex datasets.

## 5. Stages of pre-processed data

![image-20200716230850303](https://i.loli.net/2020/07/16/5b197cIqkmVWgK3.png)

<center> table 1</center>



### five stages of data processing

1. raw data
2. normalized data
3. corrected data
4. feature-selected data
5. demensionality-reduced data

### three pre-processing layers

1. measured data
2. corrected data
3. reduced data

- reduced data are used for exploratory methods that require data summaries (visualization, neighbourhood graph inference, clustering) and for computationally complex downstream analysis tools (trajectory inference)

### Pitfalls & recommendations

Use measured data for statistical testing, corrected data for visual comparison of data and reduced data for other downstream analysis based on finding the underlying biological data manifold.

# Downstream analysis

![image-20200717143720968](https://i.loli.net/2020/07/17/lFOf4H8kgnVIJ2h.png)

 <center>figure 5</center>

## Cell-level

![image-20200717164108354](https://i.loli.net/2020/07/17/GpPzoiWVDrxOKbF.png)

<center > figure 6</center>

### Cluster analysis

#### Clustering

Community detection is often faster than clustering as only neighbouring cell pairs have to be considered as belonging to the same cluster.

##### Clustering algorithm

- k-means clustering algorithm

- Based directly on a distance matrix

- Alternatives to standard Euclidean distances include cosine similarity (Haghverdi et al, 2018), correlation-based distance metrics (Kim et al, 2018) or the SIMLR method, which learns a distance metric for each dataset using Gaussian kernels (Wang et al, 2017).
- A recent comparison has suggested that **correlation-based** distances may outperform other distance metrics when used with k-means or as the basis for Gaussian kernels (Kim et al, 2018).
- Cells are assigned to clusters by minimizing intracluster distances or finding dense regions in the reduced expression space.

##### Community detection algorithm

- graph-partitioning algorithms and thus rely on a graph representation of single-cell data
- This graph representation is obtained using a **K-Nearest Neighbour approach** (KNN graph)

- multi-resolution modularity optimization (Newman & Girvan, 2004; Reichardt & Bornholdt, 2006) as implemented in the **Louvain algorithm (Blondel et al, 2008) on single-cell KNN graphs**

##### Pitfalls & recommendations

- We recommend clustering by Louvain community detection on a single-cell KNN graph.

- Clustering does not have to be performed at a single resolution. Subclustering particular cell clusters is a valid approach to focus on more detailed substructures in a dataset.

#### Cluster annotation

- it is best to use the term “cell identities” rather than “cell types”

##### two aspect to be considered

- <font color=red>Firstly, the P-values obtained for marker genes are based on the assumption that the obtained cell clusters represent the biological ground truth</font>
- Secondly, marker genes differentiate a cluster from others in the dataset and are thus dependent not only on the cell cluster, but also on the dataset composition.





##### two ways to use reference database 

###### using data-derived marker genes

- Reference webtools such as www.mousebrain.org (Zeisel et al, 2018) or http://dropviz.org/ (Saunders et al, 2018) allow users to visualize the expression of dataset marker genes in the reference dataset to facilitate cell-identity annotation.

  
###### using full gene expression profiles

- By directly comparing the gene expression profiles of annotated
  reference clusters to individual cells
- *scmap (Kiselev et al, 2018b) or Garnett (preprint: Pliner et al, 2019)*

- these methods can perform annotation and cluster assignment simultaneously,

- As cell type and state compositions differ between experimental conditions (Segerstolpe et al, 2016; Tanay & Regev, 2017), clustering based on reference data should not replace the data-driven approach.

- for large datasets that contain many clusters, the current best practice is a combination of both approaches.

   In the interest of speed, automated cell-identity annotation can be used to coarsely label cells and identify where subclustering may be needed. Subsequently, marker genes should be calculated for the dataset clusters and compared to known marker gene sets from the reference dataset or literature.

##### Pitfalls & recommendations

- Do not use marker gene P-values to validate a cell-identity cluster, especially when the detected marker genes do not help to annotate the community. P-values may be inflated.

- Note that marker genes for the same cell-identity cluster may differ between datasets purely due to dataset cell type and state compositions

- If relevant reference atlases exist, we recommend using automated cluster annotation combined with data-derived marker-gene-based manual annotation to annotate clusters.

#### Compositional analysis

- In the absence of dedicated tools, visual comparison of compositional data can be informative of changes in compositions between samples (Fig 6C).
- Future developments in this field will likely borrow from the mass cytometry (e.g. Tibshirani et al, 2002; Arvaniti & Claassen, 2017; Lun et al, 2017; Weber et al, 2018) or the microbiome literature (Gloor et al, 2017), where compositional data analysis has received more attention.

##### Pitfalls and recommendations

- Consider that statistical tests over changes in the proportion of a cell identity cluster between samples are dependent on one another.

### Trajectory analysis

![image-20200717164008984](https://i.loli.net/2020/07/17/BLDCkJ7WtUxuPl2.png)

<center> figure 7</center>

#### Trajectory inference(TI)



![image-20200718182007541](D:\OneDrive - zju.edu.cn\learning\note\Current best practices in single-cell RNA-seq analysis a tutorial.assets\image-20200718182007541.png)

<center>figure 8</center>

​                 <u>*A comparison of single-cell trajectory inference methods*</u>

- Currently available TI methods differ in the complexity of the paths that are modelled. Models range from simple linear or bifurcating trajectories, to complex graphs, trees, or multifurcating trajectories.

- no individual method performs optimally for all types of trajectories.
- Generally, any inferred trajectory should be confirmed with an alternative method to avoid method bias.
- Clusters in inferred trajectories may represent stable or metastable states (see “Metastable states”; Fig 7B and C). Subsequently, RNA velocities can be overlayed  onto the trajectory to add directionality (La Manno et al, 2018).

##### Tools

- **Slingshot** (Street et al, 2018) outperformed other methods for simple trajectories that range from linear to bi- and multifurcating models. 
- If more complex trajectories are expected, **PAGA** (Wolf et al, 2019 was recommended by the authors. 

##### Pitfalls and recommendations

- We recommend using the Saelens et al (2018)  review as a guide.

![image-20200718142151005](D:\OneDrive - zju.edu.cn\learning\note\Current best practices in single-cell RNA-seq analysis a tutorial.assets\image-20200718142151005.png)

<center> figure 9</center>

- Inferred trajectories do not have to represent a biological process. Further sources of evidence should be collected to interpret a trajectory.

#### Gene expression dynamics

##### Goal

- to support that an inferred trajectory is not the result of fitting transcriptional noise

- Potential regulatory genes are obtained by performing model selection for the genes’ dependence on pseudotime.

##### Notes

This DE test over pseudotime is confounded by the trajectory inference method in the same way that DE testing between clusters is confounded by the clustering method (see “Cluster annotation” section). Thus, P-values obtained in this set-up should not be regarded as an evaluation of significance.

##### Tools

- BEAM is a tool integrated into the Monocle TI pipeline (Qiu et al, 2017a), which allows for detection of branch-specific gene dynamics.

- users can opt for LineagePulse (https://github.com/YosefLab/LineagePulse), which considers dropout noise but is still in development, or write their own testing framework using the limma package (Ritchie et al, 2015) or standard R libraries. An example of this can be found in the online Slingshot tutorial (Street et al, 2018) and in Fig 7E. 

### Cell-level analysis unification

- By representing single-cell clusters as nodes, and trajectories between the clusters as edges, one can represent
  both the static and dynamic nature of the data. 

- partition-based graph abstraction tool (PAGA; Fig 7D; Wolf et al, 2019).

## Gene-level analysis

### Differential expression testing

#### Tools

1. DESeq2 (Love et al, 2014) and EdgeR (Robinson et al, 2010) in combination with weights estimated by ZINB-wave (Risso et al, 2018).

2. MAST(Finak et al, 2015) represents a potent alternative
   to weighted bulk DE tools.

3. While MAST has a 10-fold to 100-fold faster runtime than weighted bulk methods (Van den Berge et al, 2018), a further 10-fold speedup can be achieved using limma–voom (Law et al, 2014)

#### Pitfalls and recommendations

- DE testing should not be performed on corrected data (denoised or batch corrected, etc), but instead on measured data with technical covariates 
- <font color=red>Users should not rely on DE testing tools to correct models with confounded covariates. Model specification should be performed carefully ensuring a full-rank design matrix</font>
- We recommend using MAST or limma for DE testing.

### Gene set analysis

#### Database

- Biological process labels are stored in databases such as MSigDB (Liberzon et al, 2011)

- the Gene Ontology (Ashburner et al, 2000;The Gene Ontology Consortium, 2017), 

- or the pathway databases KEGG (Kanehisa et al, 2017) and Reactome (Fabregat et al, 2018).

(Enrichment of annotations on the gene list can be tested using a vast array of tools, which are reviewed and compared in Huang et al (2009) and Tarca et al (2013). )

- Ligand–receptor pair labels can obtained from the recent CellPhoneDB (Vento-Tormo et al, 2018) and used to interpret the highly expressed genes across clusters using statistical models (Zepp et al, 2017; Zhou et al, 2017; Cohen et al, 2018; Vento-Tormo et al, 2018).

### Gene regulatory networks

#### Goal

- Uncovering these regulatory interactions is the goal of gene regulatory network (GRN) inference methods.

#### Tools

*SCONE: Matsumoto et al, 2017;*
*PIDC: Chan et al, 2017;* 
*SCENIC: Aibar et al, 2017),*

- Inferring gene regulatory relationships is related to the detection of trajectory-associated regulatory genes. 

  several single-cell GRN inference methods use trajectories with mechanistic differential equation models (Ocone et al, 2015; Matsumoto et al, 2017).

#### Pitfalls and recommendations

<font color=red>Users should be wary of uncertainty in the inferred regulatory relationships.Modules of genes that are enriched for regulatory relationships will be more reliable than individual edges.</font>

# Analysis platforms

## Platforms

- Currently available platforms exist on the command line in R (McCarthy et al, 2017; Butler et al, 2018) or Python (Wolf et al, 2018), and as local applications (Patel, 2018; preprint: Scholz et al, 2018) or Web servers (Gardeux et al, 2017; Zhu et al, 2017) with graphical user interfaces (GUIs).

- Among command line platforms, Scater (McCarthy et al, 2017) and Seurat (Butler et al, 2018) easily interface with the large variety of analysis tools available via the R Bioconductor project (Huber et al, 2015).
- Scater has a particular strength in QC and pre-processing, 
- Seurat is arguably the most popular and comprehensive platform, which includes a large array of tools and tutorials. 
- scanpy (Wolf et al, 2018), a growing Pythonbased platform, which exhibits improved scaling to larger numbers

# Conclusion

Two avenues of development that are of particular interest due to their potential for disruption to analysis pipelines are **deep learning workflows and single-cell omic integration.**

## Deep learning

*scVI(Deep generative modeling for single-cell transcriptomics)*

## Single cell omic integration

Future single-cell platforms will have to be able to deal with different data sources such as DNA methylation (Smallwood et al, 2014), chromatin accessibility (Buenrostro et al, 2015), or protein abundance (Stoeckius et al, 2017), and include tools that integrate these modalities.



