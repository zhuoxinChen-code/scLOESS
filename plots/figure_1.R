library(ggplot2)
library(ggrepel)
library(cowplot)
library(colorspace)
library(ComplexHeatmap)

### Retrieve cell type and germ layer information
clu.cell.type <- data.frame(cluster=c(0,seq(40)),
                            germ.layer=c(
                              #1
                              "Ectoderm","Unknown","Ectoderm","Unknown","Mesoderm","Mesoderm",
                              #2
                              "Ectoderm","Ectoderm","Ectoderm","Periderm",
                              #3
                              "Ectoderm","Ectoderm","Mesoderm","Mesoderm",
                              #4
                              "Periderm","Ectoderm","Endoderm",
                              #5
                              "Unknown","Ectoderm","Ectoderm","Mesoderm",
                              #6
                              "Unknown","Mesoderm","Mesoderm","Ectoderm","Ectoderm",
                              #7
                              "Mesoderm","Ectoderm","Ectoderm","Unknown",
                              #8
                              "Ectoderm","Ectoderm","Endoderm","Ectoderm","Endoderm",
                              #9
                              "Mesoderm","Mesoderm","Ectoderm","Mesoderm","Ectoderm","Mesoderm"),
                            cell.type=c(
                              #1
                              "Neurons A (Brain)","Unknown","Neurons B","Fibroblasts A","Fibroblasts B (Somites)","Fibroblasts C \n     (Somites and vessels)",
                              #2
                              "Epidermal cells A","Neurons C (Proliferating)","Neuronal precursors","Epidermal cells B (Superficial)",
                              #3
                              "Chondrocytes A (Otic \n      capsule and neurocranium)","Retinal cells","Muscle stem cells","Skeletal muscle cells A",
                              #4
                              "Epidermal cells C \n      (Superficial)","Glial cells A","Epithelial cells A \n      (Pharyngeal pouches)",
                              #5
                              "Epidermal cells D \n      (Nose and neuromasts)","Ganglia","Radial glial cells","Endothelial cells",
                              #6
                              "Smooth muscle cells","Skeletal muscle cells B","Chondrocytes B \n      (Sclerotome)","Inner ear cells","Chondrocytes C",
                              #7
                              "Fibroblasts D (Fin folds)","Photoreceptor cells","Chondrocytes D","Hypophysis, epiphysis \n      and nose",
                              #8
                              "Epidermal cells E (Basal)","Corneal and lens cells","Intestinal cells","Osteoblast (Cranial dermal \n      bones)","Hepatocytes",
                              #9
                              "Skeletal muscle cells C","Nephric duct cells","Glial cells B","Leukocytes","Epithelial cells B (Inner ear)","Erythrocytes"),
                            gg.layer=c(
                              #1
                              "Neuroectoderm","Unknown","Neuroectoderm","Unknown","Mesoderm","Mesoderm",
                              #2
                              "Surface ectoderm","Neuroectoderm","Neuroectoderm","Periderm",
                              #3
                              "Neural crest","Neuroectoderm","Mesoderm","Mesoderm",
                              #4
                              "Periderm","Neuroectoderm","Endoderm",
                              #5
                              "Periderm + Surface ectoderm","Neuroectoderm + Otic vesicle","Neuroectoderm","Mesoderm",
                              #6
                              "Neural crest + Mesoderm","Mesoderm","Mesoderm","Otic vesicle","Neural crest",
                              #7
                              "Mesoderm","Neuroectoderm","Neural crest","Unknown",
                              #8
                              "Surface ectoderm","Neural crest + Surface ectoderm","Endoderm","Neural crest","Endoderm",
                              #9
                              "Mesoderm","Mesoderm","Neuroectoderm","Mesoderm","Otic vesicle","Mesoderm")) %>% 
  mutate_if(is.numeric, as.character) 

hue.color=c("#a5782e",
            "#caa745",
            "#e5c639",
            "#dae38d",
            "#d1ec45",   # Neural crest 5 clusters
            
            "#d76e86",
            "#e94071",
            "#a82e4c",
            "#731f2a",
            "#da9692",
            "#e43b45",
            "#b12c22",
            "#da3c71",
            "#bd0026",
            "#ed6054",   # Neuroectoderm 10 clusters
      
            "#ea3f20",
            "#e46a36",
            "#ea5f23",
            "#e28c23",   # Other ectoderm 4 clusters
            
            "#7ad2f3",
            "#589fe9",
            "#6b91f4",   # Endoderm 3 clusters
            
            "#8db56c",
            "#4bba9b",
            "#64b53f",
            "#7ae74e",
            "#b8e6b8",
            "#5be68b",
            "#6ae7df",
            "#67eabc",
            "#47b4b5",
            "#54ba77",
            "#a7e4e0",
            "#afe786",   # Mesoderm 12 clusters
           
            "#b9c1ee",
            "#9d9ce0",   # Periderm 2 clusters
            
            "#ca45e0",
            "#bc63b5",
            "#e682dd",
            "#8f3f96",
            "#b0397a"    # NA 5 clusters
            )
  
hue.fill=c("#eb5c80","#5993e5","#7dc778","#9d88da","dimgrey")

### Combine 
exp.cell.type <- exp.tsne %>% left_join(clu.cell.type,by="cluster") %>% 
  unite("label",germ.layer,cluster,sep="_",remove=F) %>% mutate(cell_id=paste0("TAG_CB_",cell,"-1"))
exp.name <- exp.cell.type %>% distinct(label,.keep_all = T) %>% 
  unite("clu.type",cluster,cell.type,sep=": ",remove=F) %>% arrange(germ.layer,gg.layer,cluster) %>%
  add_column(color=hue.color) %>% mutate_at("cluster", as.numeric) %>% arrange(germ.layer,gg.layer,cluster)

### Compute label location
tsne.cent <- exp.tsne %>% group_by(cluster) %>% select(-cell) %>% summarize_all(median) %>% 
  mutate_if(is.factor, as.character) %>% as_tibble() %>% left_join(clu.cell.type,by="cluster")

### Arrange cell types
exp.cell.type$label <- factor(exp.cell.type$label, levels=exp.name$label)
tsne.cent$germ.layer <- factor(tsne.cent$germ.layer, levels=c("Ectoderm","Endoderm","Mesoderm","Periderm","Unknown"))

### Plot figure.1A
p1.A <- exp.cell.type %>%
  ggplot(aes(x = tSNE_1, y = tSNE_2)) + 
  geom_point(aes(color=label),alpha=0.9,size=0.1) +
  scale_color_manual(values = exp.name$color,labels = exp.name$clu.type, name = "Cell type",
                     guide = FALSE
                     #guide_legend(ncol=1,override.aes = list(size=1.2))
  ) + 
  geom_label_repel(aes(fill = germ.layer,label = cluster), data = tsne.cent, nudge_x=-0.5,#point.padding=0.05,
                   color="white",fontface = 'bold',size=2,alpha=1,label.padding=0.1,force=0.01) +
  scale_fill_manual(values = hue.fill,name="Germ layer") +
  theme_bw() +
  theme(#legend.title =element_blank() ,legend.position = "bottom",legend.box = "horizontal",
    legend.text= element_text(size=5,margin = margin(0,0,0,0,"mm"),
                              #vjust=1,
                              lineheight = 1),
    legend.spacing.y = unit(0.1, 'cm'),legend.box.spacing=unit(0, 'cm'),
    legend.position = c(.8, .25),legend.justification = c("left", "top"),
    legend.key.size = unit(0.5, "mm"),
    legend.title = element_text(size = 7, color = "dimgrey"),
    plot.margin = unit(c(0.6, 0.5, 0, 0.3), "cm"),
    #axis.line = element_line(colour = "black"),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.line=element_blank(),
    axis.title.x=element_text(vjust=1,#face = "bold"
                              family = "Helvetica"),
    axis.title.y=element_text(family = "Helvetica"),
    axis.title=element_text(size=8),
    panel.background=element_rect(colour="black",size=0.5),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(name="tSNE 2") + scale_x_continuous(name="tSNE 1")

### Plot figure.1B
### Choose marker genes
chara_genes <- c("elavl3","rtn1a","lfng","mki67","her4.1","rs1a","zwi","gfap","pde6gb","tlx2","otomp","gpx2",   #neural tube
                 "col1a2","meox1","plpp3","rbp7b","myog","cdh5","tnnt3a","casq2","nkx3.2","and2","actn3a","cdh17","fcer1gl","hbae1.1",   ###mesoderm
                 "matn4","epyc","snorc","ifitm5","scinla","hpdb","col17a1a","krtt1c19e",
                 "zgc:111983","tm4sf4","tfa",
                 "cldni","anxa1b",
                 "acta2","cldnh","scg3","capsla","ogna")

### Retrieve scaled expression
char_dotplot <- DotPlot(object = cas9_larv, features = chara_genes)
mat <- char_dotplot$data %>% select(avg.exp.scaled,id,features.plot) %>% filter(!id %in% c(1,3)) %>%
  spread(id,avg.exp.scaled) %>%
  column_to_rownames("features.plot") %>% as.matrix()

### Order cell types
colord=c("0","2","7","8","11","15","19","27","37","18","24","39",
         "4","5","12","13","20","22","23","26","35","36","38","40",
         "10","25","28","33","31","6","30","16","32","34","9","14","21","17","29"
         )

p1.B <- Heatmap(t(mat),col=diverging_hcl(120, "Blue-Red")[41:120],   ### match the color of dot plots  
                cluster_rows = F, show_column_dend = F,
                show_heatmap_legend = F,
                #row_order = sort(rownames(mat)), 
                row_order = colord,
                column_order = chara_genes,
                show_column_names = F,show_row_names = F)

### Plot dot plots for different germ layers (1C-E)
chara_genes <- c("col1a2","meox1",  #4
                 "plpp3", #5
                 "rbp7b", #12
                 "myog", #13
                 "cdh5",  #20
                 "tnnt3a","casq2", #22
                 "nkx3.2", #23
                 "and2", #26
                 "actn3a", #35
                 "cdh17", #36
                 "fcer1gl", #38
                 "hbae1.1" #40
)
char_dotplot <- DotPlot(object = cas9_larv, features = chara_genes)
p1.D <- char_dotplot$data %>% filter(pct.exp>0) %>% mutate_if(is.factor, as.character) %>% 
  left_join(clu.cell.type,by=c("id"="cluster")) %>% #filter(! germ.layer=="Unknown") %>%
  mutate(id = factor(id, levels = seq(40,0))) %>%
  mutate(features.plot = factor(features.plot, levels = chara_genes)) %>% as_tibble() %>%
  filter(germ.layer=="Mesoderm") %>% unite("code",c("id","cell.type"),sep=": ",remove = F) %>%
  mutate(code = factor(code, levels = rev(filter(exp.name,germ.layer=="Mesoderm")$clu.type))) %>% 
  ggplot(aes(x=features.plot,y=code,size=pct.exp,color=avg.exp.scaled)) + geom_point() +
  scale_radius(range = c(0, 4), name="Percent of\nexpressing",
               limits=c(0,100)) + 
  scale_color_continuous_diverging(palette = "blue-red",breaks=c(-1,0,1,2,3),limits=c(-1,3),
                                   name="Average scaled\nexpression") +
  theme_bw() +
  theme(#legend.title = element_blank(),
    legend.position = "none",
    #legend.text=element_text(size=6,margin = margin(0.3,0.3,0.3,0.3,"mm")),
    axis.text = element_text(size = 5,color="black"),
    axis.text.x = element_text(angle = 45, hjust = 1,vjust=1.1),
    axis.title.x = element_blank(),#axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    axis.title.y = element_blank(),
    axis.ticks=element_line(size=0.4),
    #axis.text.y=element_text(color=rev(filter(exp.name,germ.layer=="Ectoderm")$color)),
    plot.margin = unit(c(0.6, 0.3, 0, 0.6), "cm"),
    panel.background=element_rect(colour="black",size=0.5),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_discrete(name="Cluster",position="right") + scale_x_discrete(name="Genes")


chara_genes <- c("elavl3","rtn1a",  #0
                 "lfng",   #2
                 "mki67",#7
                 "her4.1",#8
                 "rs1a",#11
                 "zwi", #15 #37
                 "gfap",#19
                 "pde6gb", #27
                 "tlx2", #18
                 "otomp",#24
                 "gpx2" #39   
                 )
char_dotplot <- DotPlot(object = cas9_larv, features = chara_genes)
p1.C <- char_dotplot$data %>% filter(pct.exp>0) %>% mutate_if(is.factor, as.character) %>% 
  left_join(clu.cell.type,by=c("id"="cluster")) %>% #filter(! germ.layer=="Unknown") %>%
  mutate(id = factor(id, levels = seq(40,0))) %>%
  mutate(features.plot = factor(features.plot, levels = chara_genes)) %>% as_tibble() %>%
  filter(grepl("Neuroectoderm",gg.layer) | grepl("Otic",gg.layer)) %>% 
  unite("code",c("id","cell.type"),sep=": ",remove = F) %>%
  mutate(code = factor(code, levels = rev(filter(exp.name,germ.layer=="Ectoderm")$clu.type))) %>% 
  ggplot(aes(x=features.plot,y=code,size=pct.exp,color=avg.exp.scaled)) + geom_point() +
  scale_radius(range = c(0, 4), name="Percentage of\nexpressing",
               limits=c(0,100)) +  
  scale_color_continuous_diverging(palette = "blue-red",breaks=c(-1,0,1,2,3),limits=c(-1,3),
                                   name="Average scaled\nexpression") +
  theme_bw() +
  theme(#legend.title = element_blank(),
  legend.position = "bottom", #legend.box = "vertical",
  legend.text=element_text(size=6,margin = margin(0,0,0,0,"mm")),
  legend.spacing = unit(0, 'cm'),legend.box.spacing=unit(0, 'cm'),
  legend.key.size = unit(0.2, "cm"),
  legend.title = element_text(size = 6, color = "dimgrey"),
  axis.text = element_text(size = 5,color="black"),
  axis.text.x = element_text(angle = 45, hjust = 1,vjust=1.1),
  axis.title.x = element_blank(),#axis.text.y=element_blank(), axis.ticks.y=element_blank(),
  axis.title.y = element_blank(),
  axis.ticks=element_line(size=0.4),
  #axis.text.y=element_text(color=rev(filter(exp.name,germ.layer=="Ectoderm")$color)),
  plot.margin = unit(c(0.6, 0.5, 0.3, 0.6), "cm"),
  panel.background=element_rect(colour="black",size=0.5),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_discrete(name="Cluster",position="right") + scale_x_discrete(name="Genes") +
  guides(size = guide_legend(label.position = "bottom",title.hjust = 0))

chara_genes <- c("matn4", #10
                 "epyc", #25
                 "snorc", #28
                 "ifitm5",#33
                 "scinla", #31
                 "hpdb",  #6
                 "col17a1a","krtt1c19e",#30
                 
                 "zgc:111983",  #16
                 "tm4sf4", #32
                 "tfa", #34
                 
                 "cldni",#9
                 "anxa1b", #14
                 
                 "acta2", #21
                 "cldnh", #17
                 #1 #3
                 "scg3","capsla" #29
                 )
char_dotplot <- DotPlot(object = cas9_larv, features = chara_genes)
p1.E <- char_dotplot$data %>% filter(pct.exp>0) %>% mutate_if(is.factor, as.character) %>% 
  left_join(clu.cell.type,by=c("id"="cluster")) %>% #filter(! germ.layer=="Unknown") %>%
  mutate(id = factor(id, levels = seq(40,0))) %>%
  mutate(features.plot = factor(features.plot, levels = chara_genes)) %>% as_tibble() %>%
  filter(!grepl("Neuroectoderm",gg.layer),!grepl("Otic",gg.layer),germ.layer!="Mesoderm") %>% 
  filter(id !=1,id !=3) %>% unite("code",c("id","cell.type"),sep=": ",remove = F) %>%
  mutate(code = factor(code, levels = rev(exp.name$clu.type))) %>% 
  ggplot(aes(x=features.plot,y=code,size=pct.exp,color=avg.exp.scaled)) + geom_point() +
  scale_radius(range = c(0, 4), name="Percentage of\nexpressing",
               limits=c(0,100)) + 
  scale_color_continuous_diverging(palette = "blue-red",breaks=c(-1,0,1,2,3),limits=c(-1,3),
                                   name="Average scaled\nexpression") +
  theme_bw() +
  theme(#legend.title = element_blank(),
    legend.position = "none",
    #legend.text=element_text(size=6,margin = margin(0.3,0.3,0.3,0.3,"mm")),
    axis.text = element_text(size = 5,color="black"),
    axis.text.x = element_text(angle = 45, hjust = 1,vjust=1.1),
    axis.title.x = element_blank(),#axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    axis.title.y = element_blank(),
    axis.ticks=element_line(size=0.4),
    #axis.text.y=element_text(color=rev(filter(exp.name,germ.layer=="Ectoderm")$color)),
    plot.margin = unit(c(0.6, 0.5, 0.1, 0.6), "cm"),
    panel.background=element_rect(colour="black",size=0.5),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_discrete(name="Cluster",position="right") + scale_x_discrete(name="Genes")

p1.BC <- plot_grid(as.ggplot(p1.B,scale = 0.95), p1.C, labels = c("B","C"),nrow=2, rel_heights = c(1, 1.5))
p1.ABC <- plot_grid(p1.A, p1.BC, labels = c("A", ""),ncol=2, rel_widths = c(1.2, 1))
p1.DE <- plot_grid(p1.D, p1.E, labels = c("D", "E"),ncol=2, rel_widths = c(1, 1.3))
p1.total <- plot_grid(p1.ABC, p1.DE, labels = c("", ""),  nrow=2, rel_heights = c(1.3, 1))
ggsave("Figure_1.pdf", height = 5.874, width = 7.25, device=cairo_pdf)
