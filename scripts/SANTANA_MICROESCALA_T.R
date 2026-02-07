#############################SANTANA GERAL###############################
#############################MICROESCALA################################

##################PREPARANDO O AMBIENTE R##############################
#setwd("C:/Users/victor.junta/Desktop/SANTANA/ANALISES/MICROESCALA")
setwd("G:/Meu Drive/BACKUP/UFLA/CEBS/Mestrado/Dissertação/Correção/GERAL/ANALISES/MICROESCALA")

library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(vegan)
 
##################TROGLOS MICROSESCALA##########################
##CARREGANDO OS DADOS##

#TABELA BASE
brutostroglo = read.csv("Geraltroglo_quadrantes_r_limpo.csv", sep=";", header = T)%>%
  replace(is.na(.), 0)

brutostroglo$ID = paste(brutostroglo$Cave, brutostroglo$Sector, sep = "_")

#TABELA DE COMUNIDADE
comunidade=brutostroglo%>%
  select(ID, Arrhopalitidae_sp1, Clivina_sp1, Endecous_sp1, Endecous_sp2, Eukoenenia_sp1, Eusarcus.cavernicola_sp1, Kinnaridae_sp1, Ochyroceratidae_sp1, Paronellidae_sp2, Pectenoniscus.santanensis_sp1, Phaneromerum_sp1, Platyartridae_sp1, Pseudochthonius_sp1, Styloniscidae._sp1, Styloniscidae._sp2)

row.names(comunidade) = comunidade$ID
comunidade=comunidade[,-1]
write.table(comunidade, "Comunidade_T.csv", row.names = T, sep = ";")

#TABELA DE COMUNIDADE BINARIA
comm.bin = decostand(comunidade, method = "pa")


#TABELA DE DADOS AMBIENTAIS
ambiente=brutostroglo%>%
  select(ID, Sector, Cave, Basin, Hidro, Distance, Subs.Div., Shelter.Div., Trophic.Div., Shelter.Av., Trophic.Av.)

row.names(ambiente) = ambiente$ID
ambiente=ambiente[,-1]

ambiente = ambiente %>%
  mutate(across(where(is.numeric), scale))

#FILTRAR VARIAVEIS AMBIENTAIS
#tiff("ColinSGeral.tiff", units="in", width=20, height=20, res=300)
#chart.Correlation(ambiente[,-c(1,2, 3, 4)], histogram = TRUE, method = c("spearman"))
#dev.off()  
#remover subs.div.

#######DBRDA COM DISTANCIA DE JACCARD#### 

#FILTRANDO QUADRANTES ZERADOS
linhas_validas = rowSums(comm.bin) > 0

comm.bin = comm.bin[linhas_validas, ]
ambiente = ambiente[linhas_validas, ]


#MATRIX DE SIMILARIDADE DE JACCARD
diss_jaccard = vegdist(comm.bin, method = "jaccard")

#DBRDA
db0 = capscale(diss_jaccard~1, data = ambiente[,-c(1, 3, 4,6)])

db = capscale(diss_jaccard ~ ., data = ambiente[,-c(1, 3, 4,6)])

summary(db)

#TESTE GLOBAL DE SIGNIFICANCIA
anova(db)

#TESTE PARA CADA VARIAVEL
anova(db, by = "term")

#PLOT BASICO
plot(db, display = "sites",type = "points")
plot(db, display = "bp", col = "red")  # vetores ambientais

#TESTANDO PELO ORDISTEP
set.seed(123)
modelo_final <- ordistep(db0, scope = formula(db), direction = "both")

#MODELO FINAL DBRDA
db_final = capscale(diss_jaccard ~ Cave + Distance, data = ambiente[,-c(1, 3, 4,6)])

summary(db_final)

set.seed(123)
anova(db_final)

set.seed(123)
anova(db_final, by = "term")


#######GRAFICO DBRDA#### 
#SCORES SITES
sites_scores = as.data.frame(db$CCA$u[, 1:2])
colnames(sites_scores) = c("CAP1", "CAP2")

#SCORES VARIAVEIS AMBIENTAIS
env_vectors = as.data.frame(db$CCA$biplot[, 1:2])
colnames(env_vectors) = c("CAP1", "CAP2")
env_vectors$Variable = rownames(env_vectors)

#ADICIONAR INFOS AMBIENTAIS NOS SCORES
sites_scores$Cave = ambiente$Cave
sites_scores$ID = rownames(sites_scores)

#COMPRIMENTO DOS VETORES
env_vectors$Length = sqrt(env_vectors$CAP1^2 + env_vectors$CAP2^2)
scale_factor = 0.7 * max(abs(range(sites_scores[, 1:2]))) / max(env_vectors$Length)
env_vectors$CAP1_scaled = env_vectors$CAP1 * scale_factor
env_vectors$CAP2_scaled = env_vectors$CAP2 * scale_factor

#CALCULAR VARIANCIA EXPLICADA
var_cap1 = round(100 * db_final$CCA$eig[1]/sum(db$CCA$eig), 1)
var_cap2 = round(100 * db_final$CCA$eig[2]/sum(db$CCA$eig), 1)
var_total = round(100 * db_final$CCA$tot.chi/db$tot.chi, 1)

# NOMES FORMATADOS
sites_scores <- sites_scores %>%
  mutate(Cave_formatted = case_when(
    Cave == "B_Verde" ~ "Cânion da Baixa Verde Cave",
    Cave == "Cedrao" ~ "Cedrão Cave",
    Cave == "Cedro" ~ "Cedro Cave", 
    Cave == "Escrevida" ~ "Pedra Escrevida Cave",
    Cave == "Labirinto" ~ "Labirinto Cave",
    Cave == "Leao" ~ "Leão Cave",
    Cave == "Padre" ~ "Padre Cave",
    Cave == "R_Bovina" ~ "Racha Bovina Cave",
    TRUE ~ as.character(Cave)
  ))

# DEFINIR CORES MANUAIS PARA CADA CAVERNA
cores_cavernas <- c(
  "Cânion da Baixa Verde Cave" = "#1f77b4",  # Azul
  "Cedrão Cave" = "#2ca02c",                 # Laranja
  "Cedro Cave" = "#2ca02c",                   # Verde
  "Pedra Escrevida Cave" = "#d62728",         # Vermelho
  "Labirinto Cave" = "#9467bd",               # Roxo
  "Leão Cave" = "#7c562b",                    # Marrom
  "Padre Cave" = "#ff7f0e",                   # Rosa
  "Racha Bovina Cave" = "#1f77b4"             # Cinza
)

# DEFINIR POSIÇÕES MANUAIS PARA OS NOMES DAS CAVERNAS
p# DEFINIR POSIÇÕES MANUAIS PARA OS NOMES DAS CAVERNAS
posicoes_nomes <- data.frame(
  Cave_formatted = c(
    "Cânion da Baixa Verde Cave",
    "Cedrão Cave", 
    "Cedro Cave",
    "Pedra Escrevida Cave",
    "Labirinto Cave",
    "Leão Cave",
    "Padre Cave",
    "Racha Bovina Cave"
  ),
  x_pos = c(0.34, -0.205, -0.185, 0.28, 0.209, -0.03, 0.01, 0.33),
  y_pos = c(0.30, 0.14, 0.14, 0.13, -0.44, -0.09, 0.05, 0.37)
) %>%
  dplyr::mutate(ID = case_when(
    Cave_formatted == "Cânion da Baixa Verde Cave" ~ 1,
    Cave_formatted == "Cedrão Cave"                ~ 16,
    Cave_formatted == "Cedro Cave"                 ~ 15,
    Cave_formatted == "Pedra Escrevida Cave"       ~ 5,
    Cave_formatted == "Labirinto Cave"             ~ 3,
    Cave_formatted == "Leão Cave"                  ~ 22,
    Cave_formatted == "Padre Cave"                 ~ 2,
    Cave_formatted == "Racha Bovina Cave"          ~ 10,
    TRUE ~ NA_integer_
  ))

# GERANDO O GRÁFICO
ggplot() +
  
  # PLOTAR ELIPSES PARA CADA CAVERNA
  stat_ellipse(data = sites_scores, 
             aes(x = CAP1, y = CAP2, color = Cave_formatted),
             type = "norm", level = 0.95, linewidth = 0.8, alpha = 0.6) +
  
  # PLOTAR OS PONTOS DOS QUADRANTES
  geom_point(data = sites_scores, 
             aes(x = CAP1, y = CAP2, color = Cave_formatted),
             size = 6, alpha = 0.7) +
  
  # PLOTAR VETORES
  geom_segment(data = env_vectors[env_vectors$Variable == "Distance",],
               aes(x = 0, y = 0, xend = CAP1_scaled, yend = CAP2_scaled),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "black", linewidth = 1.2) +
  
  # ROTULAR VETORES
  geom_text_repel(data = env_vectors[env_vectors$Variable == "Distance",],
                  aes(x = CAP1_scaled * 1.1, y = CAP2_scaled * 1.1, 
                      label = Variable),
                  color = "black", fontface = "bold", size = 6,
                  box.padding = 0.4, segment.color = "black", segment.alpha = 0.3) +
  
  # ADICIONAR NOMES DAS CAVERNAS COM POSIÇÕES MANUAIS
  geom_text(data = posicoes_nomes,
            aes(x = x_pos, y = y_pos, label = ID, color = Cave_formatted),
            size = 6, alpha = 1, fontface = "bold", check_overlap = TRUE) +
  
  # LINHAS DE REFERENCIA
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
  
  # FORMATACAO DO GRAFICO
  labs(title = "Microscale - dbRDA - Troglobitic species",
       subtitle = paste("Jaccard distance | ", var_total, "% explained"),
       x = paste("CAP1 (", var_cap1, "%)"),
       y = paste("CAP2 (", var_cap2, "%)")) +
  
  # ESCALA DE CORES MANUAL
  scale_color_manual(values = cores_cavernas) +
  
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.9),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(face = "bold", size = 16),
        legend.position = "none")

ggsave("DBRDA_Troglos_Quad.tiff", width = 14, height = 10, dpi = 300, bg = "white")

###################################################################################
#########################ENTENDENDO AS RELAÇÕES ESPACIAIS#####################

###COORDENADAS###
coord_quad = read.csv("quadrantes_utm_ID.csv", sep=";", header = T)
coord_quad = coord_quad[linhas_validas, ]
###MATRIX DE COORDENADAS###
coords_matrix = cbind(coord_quad$x, coord_quad$y)

dist_geo = dist(coords_matrix)

# VERIFICAR ESTRUTURA DA MATRIZ DE DISTANCIAS
cat("Dimensão da matriz de distância:", dim(as.matrix(dist_geo)), "\n")
cat("Distância mínima entre quadrantes:", min(dist_geo), "m\n")
cat("Distância máxima entre quadrantes:", max(dist_geo)/1000, "km\n")


#TESTE DE MANTEL
mantel_result = mantel(diss_jaccard, dist_geo, 
                        method = "spearman",
                        permutations = 9999)

#CONVERTER VETORES PARA PLOTAGEM
distancias = as.vector(as.matrix(dist_geo))
dissimilaridades = as.vector(as.matrix(diss_jaccard))

# Criar dataframe para plot
dados_plot = data.frame(
  Distancia_m = distancias,
  Distancia_km = distancias/1000,
  Dissimilaridade = dissimilaridades
)

# Plot da relação
ggplot(dados_plot, aes(x = Distancia_km, y = Dissimilaridade)) +
  geom_point(alpha = 0.2, size = 2, color = "gray40") +
  geom_smooth(method = "gam", color = "blue", se = TRUE, linewidth = 1.2) +
  labs(
    title = "Geographic Distance vs. Dissimilarity",
    subtitle = paste("Mantel r =", round(mantel_result$statistic, 4),
                     "| p =", mantel_result$signif),
    x = "Distance between quadrats (km)",
    y = "Jaccard dissimilarity"
  ) +
  theme_minimal(base_size = 16) +  # aumenta o tamanho de fonte geral
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  # título centralizado e maior
    plot.subtitle = element_text(hjust = 0.5, size = 14),              # subtítulo centralizado
    axis.title = element_text(face = "bold", size = 16),               # títulos dos eixos em negrito
    axis.text = element_text(size = 14),                               # rótulos dos eixos maiores
    panel.grid.minor = element_blank(),                                # remove grade menor
    panel.grid.major = element_line(color = "gray85")                  # grade principal mais suave
  )


ggsave("Mantel_troglo_Quad.tiff", width = 14, height = 10, dpi = 300, bg = "white")


######################PARTIÇÃO DA VARIAÇÃO###################################
###DADOS AMBIENTAIS###
X_env = ambiente[, c("Distance")]

###GERANDO POLINOMIOS DE TENDENCIA ESPACIAL###

#COORDENADAS
coords = coord_quad[, c("x", "y")]
rownames(coords) = coord_quad$ID

#CALCULAR MATRIZ DE DISTANCIAS EUCLIDIANA GEOGRAFICA
dist_geo = dist(coords)

#GERAR AUTOVETORES 
library(adespatial)
dbmem = dbmem(dist_geo, silent = FALSE) 

#VISUALIZAR AUTOVALORES
plot(attributes(dbmem)$values,
     xlab = "Ordem do dbMEM", ylab = "Autovalor",
     main = "Quebra-de-espelho para seleção de dbMEMs")
abline(h = 0, col = "red")

#CRIAR DATAFRAME DOS AUTOVALORES
X_spa = as.data.frame(dbmem)

#PARTIÇÃO DA VARIANCIA
vp = varpart(comm.bin, ~ Cave + Distance, # Fórmula para Environment
             X_spa, # Matriz para Space
             data = ambiente)

#VISUALIZAR
vp
plot(vp,
     digits = 2,
     bg = c("skyblue", "pink"),
     Xnames = c("Ambiente", "Espaço"))

#SIGNIFICANCIA DA FRAÇÃO PURAMENTE AMBIENTAL
rda_env = rda(comm.bin ~ Cave + Distance + Condition(MEM1 + MEM2),
              data = cbind(ambiente, X_spa))

set.seed(123)
anova(rda_env, by = "margin")

set.seed(123)
anova(rda_env)

#SIGNIFICANCIA DA FRAÇÃO PURAMENTE ESPACIAL
rda_spa = rda(comm.bin ~ MEM1 + MEM2  + Condition(Cave + Distance),
              data = cbind(ambiente, X_spa))
set.seed(123)
anova(rda_spa) 

set.seed(125)
anova(rda_spa, by = "margin") 


#DIAGRAMA DE VENN
library(eulerr) 

#PREPARAR OS DADOS
r2_env = vp$part$fract$R.squared[1]    # [a+c] - Variação explicada pelo Ambiente
r2_spa = vp$part$fract$R.squared[2]    # [b+c] - Variação explicada pelo Espaço
r2_total = vp$part$fract$R.squared[3]  # [a+b+c] - Variação total explicada

r2_shared = r2_env + r2_spa - r2_total

r2_shared = max(0, r2_env + r2_spa - r2_total)
exclusive_env = max(0, r2_env - r2_shared)
exclusive_spa = max(0, r2_spa - r2_shared)

#VETORES COM AS ÁREAS PARA O VENN
fit = euler(c("Environmental" = r2_env - r2_shared,
              "Spatial" = r2_spa - r2_shared,
              "Environmental&Spatial" = r2_shared))

#PLOTAR DIAGRAMA

tiff("Venn_Micro_T.tiff", 
     width = 16,         # Largura em cm (ótimo para dissertações)
     height = 12,        # Altura em cm
     units = "cm",       # Unidade = centímetros
     res = 600,          # Alta resolução (600 dpi)
     compression = "lzw") # Compressão sem perda de qualidade


plot(fit,
     quantities = list(type = "percent", fontsize = 10),
     fills = list(fill = c("#1f77b4", "#ff7f0e", "#2ca02c"), alpha = 0.7),
     labels = list(font = 2),
     main = "Variance Partitioning")

dev.off()


###########################################################################################################################
#####################################ANÁLISE DE RIQUEZA###################################################################
##########################################################################################################################

##################PREPARANDO O AMBIENTE R##############################
setwd("C:/Users/victor.junta/Desktop/SANTANA/ANALISES/MICROESCALA")
#setwd("G:/Meu Drive/BACKUP/UFLA/CEBS/Mestrado/Dissertação/Correção/GERAL/ANALISES/MICROESCALA")

library(dplyr)
library(ggplot2)
library(lme4)        # para GLMM
library(lmerTest)    # p-valor nos modelos mistos
library(MASS)        # glm.nb
library(car)         # teste de Anova
library(performance) # diagnóstico do modelo
library(vegan)
library(stringr)


#TABELA BASE
brutostroglo = read.csv("Geraltroglo_quadrantes_r_limpo.csv", sep=";", header = T)%>%
  replace(is.na(.), 0)

brutostroglo$ID = paste(brutostroglo$Cave, brutostroglo$Sector, sep = "_")

#TABELA DE COMUNIDADE
comunidade=brutostroglo%>%
  dplyr::select(ID, Arrhopalitidae_sp1 : Styloniscidae._sp2)

row.names(comunidade) = comunidade$ID
comunidade=comunidade[,-1]

#TABELA DE DADOS AMBIENTAIS
ambiente=brutostroglo%>%
  dplyr::select(ID, Sector, Cave, Basin, Hidro, Distance, Subs.Div., Shelter.Div., Trophic.Div., Shelter.Av., Trophic.Av.)

row.names(ambiente) = ambiente$ID


ambiente = ambiente %>%
  mutate(across(where(is.numeric), scale))

ambiente$Sector = str_remove(ambiente$ID, "(?i)Q\\d+.*$") 
####CALCULANDO A RIQUEZA DOS QUADRANTES####
riqueza = (comunidade > 0) %>%  
  rowSums() %>%                 
  as.data.frame()

colnames(riqueza) = "Riqueza"
riqueza$ID = rownames(comunidade)



####JUNTANDO COM OS DADOS AMBIENTAIS####
dados_riqueza = left_join(riqueza, ambiente, by = "ID")

####EXPLORAÇÃO INICIAL####
summary(dados_riqueza$Riqueza)

#HISTOGRAMA
hist(dados_riqueza$Riqueza, main = "Distribuição da riqueza", xlab = "Número de espécies")

#SHAPIRO TEST PARA NORMALIDADE
shapiro.test(dados_riqueza$Riqueza) 

####CORRELAÇÃO ESPACIAL DA RIQUEZA COM MANTEL####

#DISTANCIA EUCLIDIANA DA RIQUEZA
riqueza_dist = dist(dados_riqueza$Riqueza) 

###COORDENADAS###
coord_quad = read.csv("quadrantes_utm_ID.csv", sep=";", header = T)

###MATRIX DE COORDENADAS###
coords_matrix = cbind(coord_quad$x, coord_quad$y)

dist_geo = dist(coords_matrix)

###TESTE DE MANTEL###
set.seed(123)
mantel_riqueza = mantel(riqueza_dist, dist_geo, method = "spearman", permutations = 9999)
mantel_riqueza

###MODELOS
#MODELO GLMM
m1 = glmer(Riqueza ~ Distance + Shelter.Div.+Trophic.Div.+Shelter.Av.+Trophic.Av.+(1|Sector:Cave), data = dados_riqueza, family = "poisson")

summary(m1)

#VERIFICAR OVERDISPERSION
check_overdispersion(m1)

#TESTAR SIGNIFICANCIA
Anova(m1, type = "II")

###VISUALIZAR OS EFEITOS
library(ggeffects)

#EFEITO DA DISTANCIA
pred_dist = ggpredict(m1, terms = "Distance")
plot(pred_dist) + labs(x = "Distance", y = "Predicted richness")

#SALVANDO PREDICT
pred_troglo <- pred_dist  
saveRDS(pred_troglo, "pred_troglo.rds")

#SALVANDO DADOS
dados_t=dados_riqueza%>%
  dplyr::select(Riqueza,Distance)
dados_t$Grupo <- "Troglo"
saveRDS(dados_t, "dados_t.rds")

