#############################SANTANA GERAL###############################
#############################MESOESCALA################################

##################PREPARANDO O AMBIENTE R##############################
#setwd("C:/Users/victor.junta/Desktop/SANTANA/ANALISES/MESOESCALA")
setwd("G:/Meu Drive/BACKUP/UFLA/CEBS/Mestrado/Dissertação/Correção/GERAL/ANALISES/MESOESCALA/N_TROGLO")

library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(vegan)
library(PerformanceAnalytics)

##################NÃO TROGLOS MICROESCALA##########################
##CARREGANDO OS DADOS##

#TABELA BASE
brutosntroglo = read.csv("Geral_ntroglo_setores_r_limpo.csv", sep=";", header = T) %>%
  replace(is.na(.), 0)

brutosntroglo$ID = paste(brutosntroglo$Cave, brutosntroglo$Sector, sep = "_")

#TABELA DE COMUNIDADE
comunidade_nt = brutosntroglo %>%
  select(ID, Blattodea_sp4 : Novamundoniscus_sp1)

row.names(comunidade_nt) = comunidade_nt$ID
comunidade_nt = comunidade_nt[,-1]
write.table(comunidade_nt, "Comunidade_NT_Mesoescala.csv", row.names = T, sep = ";")

#TABELA DE COMUNIDADE BINARIA
comm.bin_nt = decostand(comunidade_nt, method = "pa")

#TABELA DE DADOS AMBIENTAIS
ambiente_nt = brutosntroglo %>%
  select(ID, Sector, Cave, Humidity, Temperature, Distance, Subs.Div., 
         Shelter.Div., Trophic.Div., Shelter.Av., Trophic.Av.)

row.names(ambiente_nt) = ambiente_nt$ID
ambiente_nt = ambiente_nt[,-1]

ambiente_nt <- ambiente_nt %>%
  mutate(across(where(is.numeric), scale))

#FILTRAR VARIAVEIS AMBIENTAIS
tiff("ColinSGeral.tiff", units="in", width=20, height=20, res=300)
chart.Correlation(ambiente[,-c(1:5)], histogram = TRUE, method = c("spearman"))
dev.off()  
#remover subs.div.
#######DBRDA COM DISTANCIA DE JACCARD#### 

#FILTRANDO QUADRANTES ZERADOS
linhas_validas_nt = rowSums(comm.bin_nt) > 0

comm.bin_nt = comm.bin_nt[linhas_validas_nt, ]
ambiente_nt = ambiente_nt[linhas_validas_nt, ]

#MATRIX DE SIMILARIDADE DE JACCARD
diss_jaccard_nt = vegdist(comm.bin_nt, method = "jaccard")

#DBRDA
db0_nt = capscale(diss_jaccard_nt ~ 1, data = ambiente_nt[,-c(1)])
db_nt = capscale(diss_jaccard_nt ~ ., data = ambiente_nt[,-c(1)])

summary(db_nt)

#TESTES
set.seed(123)
anova(db_nt)

set.seed(123)
anova(db_nt, by = "term")

#ORDISTEP
set.seed(123)
modelo_final_nt <- ordistep(db0_nt, scope = formula(db_nt), direction = "both")

#MODELO FINAL DBRDA 
db_final_nt = capscale(diss_jaccard_nt ~ Cave + Distance + Shelter.Av. + Humidity, data = ambiente_nt[,-c(1)])

summary(db_final_nt)

set.seed(123)
anova(db_final_nt)

set.seed(123)
anova(db_final_nt, by="term")

#######GRAFICO DBRDA#### 
sites_scores_nt = as.data.frame(db_final_nt$CCA$u[, 1:2])
colnames(sites_scores_nt) = c("CAP1", "CAP2")

env_vectors_nt = as.data.frame(db_final_nt$CCA$biplot[, 1:2])
colnames(env_vectors_nt) = c("CAP1", "CAP2")
env_vectors_nt$Variable = rownames(env_vectors_nt)

sites_scores_nt$Cave = ambiente_nt$Cave
sites_scores_nt$ID = rownames(sites_scores_nt)

# Escalando vetores
env_vectors_nt$Length = sqrt(env_vectors_nt$CAP1^2 + env_vectors_nt$CAP2^2)
scale_factor_nt = 1.2 * max(abs(range(sites_scores_nt[, 1:2]))) / max(env_vectors_nt$Length)
env_vectors_nt$CAP1_scaled = env_vectors_nt$CAP1 * scale_factor_nt
env_vectors_nt$CAP2_scaled = env_vectors_nt$CAP2 * scale_factor_nt

# Variância explicada
var_cap1_nt = round(100 * db_final_nt$CCA$eig[1]/sum(db_nt$CCA$eig), 1)
var_cap2_nt = round(100 * db_final_nt$CCA$eig[2]/sum(db_nt$CCA$eig), 1)
var_total_nt = round(100 * db_final_nt$CCA$tot.chi/db_nt$tot.chi, 1)

#nomes nonitos
sites_scores_nt <- sites_scores_nt %>%
  mutate(Cave_formatted = case_when(
    Cave == "B_Verde"     ~ "Cânion da Baixa Verde Cave",
    Cave == "Padre"       ~ "Padre Cave",
    Cave == "Labirinto"   ~ "Labirinto Cave",
    Cave == "Boqueirao"   ~ "Boqueirão Cave",
    Cave == "Escrevida"   ~ "Pedra Escrevida Cave",
    Cave == "2Cobras"     ~ "Duas Cobras Cave",
    Cave == "Tunel2"      ~ "Tunel 2 Cave",
    Cave == "S_Geraldo"   ~ "São Geraldo Cave",
    Cave == "Cumbra"      ~ "Olho D'Água do Cumbra Cave",
    Cave == "R_Bovina"    ~ "Racha Bovina Cave",
    Cave == "Tunel1"      ~ "Tunel 1 Cave",
    Cave == "C_Flor"      ~ "Couve Flor Cave",
    Cave == "G_Cruz"      ~ "Geraldo Cruz Cave",
    Cave == "F_Obliqua"   ~ "Fenda Oblíquoa Cave",
    Cave == "Cedro"       ~ "Cedro Cave",
    Cave == "Cedrao"      ~ "Cedrão Cave",
    Cave == "Cedriculo"   ~ "Cedrículo Cave",
    Cave == "Pajeu"       ~ "Pajeú Cave",
    Cave == "Cristal"     ~ "Cristal Cave",
    Cave == "Salobro"     ~ "Salobro Cave",
    Cave == "Grota"       ~ "Grota Cave",
    Cave == "Leao"        ~ "Leão Cave",
    Cave == "Cinquentona" ~ "Cinquentona Cave",
    Cave == "Escrevidinha"~ "Escrevidinha Cave",
    TRUE ~ as.character(Cave)
  ),
  ID = case_when(
    Cave_formatted == "Cânion da Baixa Verde Cave" ~ 1,
    Cave_formatted == "Padre Cave"                 ~ 2,
    Cave_formatted == "Labirinto Cave"             ~ 3,
    Cave_formatted == "Boqueirao Cave"             ~ 4,
    Cave_formatted == "Pedra Escrevida I Cave"     ~ 5,
    Cave_formatted == "Duas Cobras Cave"           ~ 6,
    Cave_formatted == "Tunel 2 Cave"               ~ 7,
    Cave_formatted == "São Geraldo Cave"           ~ 8,
    Cave_formatted == "Olho D'Água do Cumbra Cave"  ~ 9,
    Cave_formatted == "Racha Bovina Cave"          ~ 10,
    Cave_formatted == "Tunel 1 Cave"               ~ 11,
    Cave_formatted == "Couve Flor Cave"            ~ 12,
    Cave_formatted == "Geraldo Cruz Cave"          ~ 13,
    Cave_formatted == "Fenda Obliquoa Cave"         ~ 14,
    Cave_formatted == "Cedro Cave"                 ~ 15,
    Cave_formatted == "Cedrão Cave"                ~ 16,
    Cave_formatted == "Cedrículo Cave"             ~ 17,
    Cave_formatted == "Pajeú Cave"                 ~ 18,
    Cave_formatted == "Cristal Cave"               ~ 19,
    Cave_formatted == "Salobro Cave"               ~ 20,
    Cave_formatted == "Grota Cave"                 ~ 21,
    Cave_formatted == "Leão Cave"                  ~ 22,
    Cave_formatted == "Cinquentona Cave"           ~ 23,
    Cave_formatted == "Escrevidinha Cave"          ~ 24,
    TRUE ~ NA_integer_
  ))

# DEFINIR CORES MANUAIS PARA CADA CAVERNA
cores_cavernas <- c(
  "Cânion da Baixa Verde Cave" = "#9467bd",   # Azul
  "Cedrão Cave"                = "#d62728",   # Verde (mantida da original)
  "Cedro Cave"                 = "#d62728",   # Verde (mantida da original)
  "Pedra Escrevida Cave"       = "#d62728",   # Vermelho
  "Labirinto Cave"             = "#9467bd",   # Roxo
  "Leão Cave"                  = "#d62728",   # Marrom (mantida da original)
  "Padre Cave"                 = "#ff7f0e",   # Laranja (mantida da original)
  "Racha Bovina Cave"          = "#9467bd",   # Azul (mantida da original)
  "Boqueirão Cave"             = "#9467bd",   # Verde-amarelado
  "Couve Flor Cave"            = "#9467bd",   # Azul-petróleo
  "São Geraldo Cave"           = "#d62728",   # Laranja-claro
  "Tunel 2 Cave"               = "#d62728",   # Verde-claro
  "Geraldo Cruz Cave"          = "#9467bd",   # Vermelho-claro
  "Duas Cobras Cave"           = "#9467bd",   # Lilás
  "Olho D'Água do Cumbra Cave" = "#9467bd",   # Bege
  "Fenda Oblíquoa Cave"        = "#d62728",   # Rosa-claro
  "Tunel 1 Cave"               = "#d62728",   # Amarelo-claro
  "Pajeú Cave"                 = "#9467bd",   # Azul-claro
  "Grota Cave"                 = "#ff7f0e",   # Azul-escuro
  "Salobro Cave"               = "#ff7f0e",   # Verde-escuro
  "Cedrículo Cave"             = "#d62728",   # Marrom-claro
  "Cinquentona Cave"           = "#d62728",   # Vinho
  "Escrevidinha Cave"          = "#d62728",   # Roxo-escuro
  "Cristal Cave"               = "#d62728"    # Azul-claro pastel
)


# DEFINIR POSIÇÕES MANUAIS PARA OS NOMES DAS CAVERNAS
posicoes_nomes <- data.frame(
  Cave_formatted = c(
    "Boqueirão Cave",
    "Cânion da Baixa Verde Cave",
    "Cedrão Cave",
    "Cedro Cave",
    "Cedrículo Cave",
    "Cinquentona Cave",
    "Cristal Cave",
    "Couve Flor Cave",
    "Duas Cobras Cave",
    "Escrevidinha Cave",
    "Fenda Oblíquoa Cave",
    "Geraldo Cruz Cave",
    "Grota Cave",
    "Labirinto Cave",
    "Leão Cave",
    "Olho D'Água do Cumbra Cave",
    "Padre Cave",
    "Pajeú Cave",
    "Pedra Escrevida Cave",
    "Racha Bovina Cave",
    "Salobro Cave",
    "São Geraldo Cave",
    "Tunel 1 Cave",
    "Tunel 2 Cave"
  ),
  x_pos = c(0.04702866, 0.07890073, 0.06494439, -0.00295923, 0.1024315,
            0.18748468,  0.11508384,  0.07890073,  0.08404293, 0.14580310,
            0.06408638,  0.0621592, -0.04659518,  0.14605454,  0.07002055,
            0.06106332, -0.10589400,  0.04469512,  0.01858766,  0.13011496,
            -0.08275971,  0.11146632,  0.14611016,  0.07022949),
  y_pos = c(0.147904912, 0.283276539, -0.184898257, -0.136140908, -0.184898257,
            -0.077412522, -0.078296280,  0.023276539, 0.075384970, -0.058236328,
            -0.149556250,  0.219896143,  0.077273478, 0.020770866, -0.061173233,
            0.168674336, -0.048988053,  0.090427900, -0.051497042, -0.022507128,
            0.081422726, -0.130977245, -0.171645315, -0.114741590)
) %>%
  dplyr::mutate(ID = case_when(
    Cave_formatted == "Cânion da Baixa Verde Cave" ~ 1,
    Cave_formatted == "Padre Cave"                 ~ 2,
    Cave_formatted == "Labirinto Cave"             ~ 3,
    Cave_formatted == "Boqueirão Cave"             ~ 4,
    Cave_formatted == "Pedra Escrevida Cave"       ~ 5,
    Cave_formatted == "Duas Cobras Cave"           ~ 6,
    Cave_formatted == "Tunel 2 Cave"               ~ 7,
    Cave_formatted == "São Geraldo Cave"           ~ 8,
    Cave_formatted == "Olho D'Água do Cumbra Cave" ~ 9,
    Cave_formatted == "Racha Bovina Cave"          ~ 10,
    Cave_formatted == "Tunel 1 Cave"               ~ 11,
    Cave_formatted == "Couve Flor Cave"            ~ 12,
    Cave_formatted == "Geraldo Cruz Cave"          ~ 13,
    Cave_formatted == "Fenda Oblíquoa Cave"        ~ 14,
    Cave_formatted == "Cedro Cave"                 ~ 15,
    Cave_formatted == "Cedrão Cave"                ~ 16,
    Cave_formatted == "Cedrículo Cave"             ~ 17,
    Cave_formatted == "Pajeú Cave"                 ~ 18,
    Cave_formatted == "Cristal Cave"               ~ 19,
    Cave_formatted == "Salobro Cave"               ~ 20,
    Cave_formatted == "Grota Cave"                 ~ 21,
    Cave_formatted == "Leão Cave"                  ~ 22,
    Cave_formatted == "Cinquentona Cave"           ~ 23,
    Cave_formatted == "Escrevidinha Cave"          ~ 24,
    TRUE ~ NA_integer_
  ))


# Gráfico
ggplot() +
  stat_ellipse(data = sites_scores_nt, 
               aes(x = CAP1, y = CAP2, color = Cave_formatted),
               type = "norm", level = 0.95, linewidth = 0.8, alpha = 0.6) +
  geom_point(data = sites_scores_nt, 
             aes(x = CAP1, y = CAP2, color = Cave_formatted),
             size = 6, alpha = 0.7) +
  geom_segment(data = env_vectors_nt[env_vectors_nt$Variable == "Humidity",],
               aes(x = 0, y = 0, xend = CAP1_scaled*0.9, yend = CAP2_scaled*0.9),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "black", linewidth = 1.2) +
  geom_segment(data = env_vectors_nt[env_vectors_nt$Variable == "Distance",],
               aes(x = 0, y = 0, xend = CAP1_scaled, yend = CAP2_scaled),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "black", linewidth = 1.2) +
  geom_segment(data = env_vectors_nt[env_vectors_nt$Variable == "Shelter.Av.",],
               aes(x = 0, y = 0, xend = CAP1_scaled, yend = CAP2_scaled),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "black", linewidth = 1.2) +
  geom_text_repel(data = env_vectors_nt[env_vectors_nt$Variable == "Humidity",],
                  aes(x = CAP1_scaled * 1, y = CAP2_scaled * 1, 
                      label = Variable),
                  color = "black", fontface = "bold", size = 6) +
  geom_text_repel(data = env_vectors_nt[env_vectors_nt$Variable == "Distance",],
                  aes(x = CAP1_scaled * 1.1, y = CAP2_scaled * 1.1, 
                      label = Variable),
                  color = "black", fontface = "bold", size = 6) +
  geom_text_repel(data = env_vectors_nt[env_vectors_nt$Variable == "Shelter.Av.",],
                  aes(x = CAP1_scaled * 1.2, y = CAP2_scaled * 1.3, 
                      label = Variable),
                  color = "black", fontface = "bold", size = 6) +
  geom_text(data = posicoes_nomes,
            aes(x = x_pos, y = y_pos, label = ID, color = Cave_formatted),
            size = 6, alpha = 1, fontface = "bold", check_overlap = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
  labs(title = "Mesoscale - dbRDA - Non-troglobitic species",
       subtitle = paste("Jaccard distance | ", var_total_nt, "% explained"),
       x = paste("CAP1 (", var_cap1_nt, "%)"),
       y = paste("CAP2 (", var_cap2_nt, "%)")) +
  scale_color_manual(values = cores_cavernas) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.9),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(face = "bold", size = 16),
        legend.position = "none")


ggsave("DBRDA_NaoTroglos_Meso.tiff", width = 14, height = 10, dpi = 300, bg = "white")

###################################################################################
#########################ENTENDENDO AS RELAÇÕES ESPACIAIS#####################
brutosntroglo$ID


coord_setor_nt = read.csv("setores_utm_ID.csv", sep=";", header = T)
coord_setor_nt <- coord_setor_nt[coord_setor_nt$ID %in% rownames(ambiente_nt), ]
coord_setor_nt$ID

coords_matrix_nt = cbind(coord_setor_nt$x, coord_setor_nt$y)
dist_geo_nt = dist(coords_matrix_nt)

mantel_result_nt = mantel(diss_jaccard_nt, dist_geo_nt, 
                          method = "spearman",
                          permutations = 9999)

distancias_nt = as.vector(as.matrix(dist_geo_nt))
dissimilaridades_nt = as.vector(as.matrix(diss_jaccard_nt))

dados_plot_nt = data.frame(
  Distancia_m = distancias_nt,
  Distancia_km = distancias_nt/1000,
  Dissimilaridade = dissimilaridades_nt
)

ggplot(dados_plot_nt, aes(x = Distancia_km, y = Dissimilaridade)) +
  geom_point(alpha = 0.2, size = 2, color = "gray40") +
  geom_smooth(method = "gam", color = "blue", se = TRUE, linewidth = 1.2) +
  labs(
    title = "Geographic Distance vs. Dissimilarity",
    subtitle = paste("Mantel r =", round(mantel_result_nt$statistic, 4),
                     "| p =", mantel_result_nt$signif),
    x = "Distance between transects (km)",
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


ggsave("Mantel_NaoTroglos_Setor.tiff", width = 14, height = 10, dpi = 300, bg = "white")


######################PARTIÇÃO DA VARIAÇÃO###################################
###DADOS AMBIENTAIS###
X_env = ambiente_nt[, c("Distance", "Shelter.Av.", "Humidity")]

###GERANDO POLINOMIOS DE TENDENCIA ESPACIAL###

#COORDENADAS
coords = coord_setor_nt[, c("x", "y")]
rownames(coords) = coord_setor_nt$ID

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
vp = varpart(comm.bin_nt, ~ Cave + Distance + Shelter.Av. + Humidity, # Fórmula para Environment
             X_spa, # Matriz para Space
             data = ambiente_nt)

#VISUALIZAR
vp
plot(vp,
     digits = 2,
     bg = c("skyblue", "pink"),
     Xnames = c("Ambiente", "Espaço"))

#SIGNIFICANCIA DA FRAÇÃO PURAMENTE AMBIENTAL
rda_env = rda(comm.bin_nt ~ Cave + Distance + Shelter.Av. + Humidity + Condition(MEM1 + MEM2),
              data = cbind(ambiente_nt, X_spa))

set.seed(125)
anova(rda_env, by = "margin")

set.seed(125)
anova(rda_env)

#SIGNIFICANCIA DA FRAÇÃO PURAMENTE ESPACIAL
rda_spa = rda(comm.bin_nt ~ MEM1 + MEM2  + Condition(Cave + Distance + Shelter.Av. + Humidity),
              data = cbind(ambiente_nt, X_spa))
set.seed(123)
anova(rda_spa) 

set.seed(123)
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

tiff("Venn_Meso_NT.tiff", 
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
setwd("C:/Users/victor.junta/Desktop/SANTANA/ANALISES/MESOESCALA")
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
brutosntroglo = read.csv("Geral_ntroglo_setores_r_limpo.csv", sep=";", header = T)%>%
  replace(is.na(.), 0)

brutosntroglo$ID = paste(brutosntroglo$Cave, brutosntroglo$Sector, sep = "_")

#TABELA DE COMUNIDADE
comunidade=brutosntroglo%>%
  dplyr::select(ID, Blattodea_sp4 : Novamundoniscus_sp1)

row.names(comunidade) = comunidade$ID
comunidade=comunidade[,-1]

#TABELA DE DADOS AMBIENTAIS
ambiente=brutosntroglo%>%
  dplyr::select(ID, Sector, Cave, Basin, Hidro, Distance, Subs.Div., Shelter.Div., Trophic.Div., Shelter.Av., Trophic.Av., Temperature, Umidity)

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
coord_quad = read.csv("setores_utm_ID.csv", sep=";", header = T)
coord_quad <- coord_quad[coord_quad$ID %in% rownames(ambiente), ]

###MATRIX DE COORDENADAS###
coords_matrix = cbind(coord_quad$x, coord_quad$y)

dist_geo = dist(coords_matrix)

###TESTE DE MANTEL###
set.seed(123)
mantel_riqueza = mantel(riqueza_dist, dist_geo, method = "spearman", permutations = 9999)
mantel_riqueza

###MODELOS
#MODELO GLMM
m1 = glmer(Riqueza ~ Distance + Temperature + Umidity + Shelter.Div.+Trophic.Div.+Shelter.Av.+Trophic.Av.+(1|Cave), data = dados_riqueza, family = "poisson")

summary(m1)


#VERIFICAR OVERDISPERSION
check_overdispersion(m1)

#TESTAR SIGNIFICANCIA
Anova(m1, type = "II")

###VISUALIZAR OS EFEITOS
library(ggeffects)

#EFEITO DA DISTANCIA
pred_dist = ggpredict(m1, terms = c("Umidity"))
plot(pred_dist) + labs(x = "Distance", y = "Predicted richness")

#SALVANDO PREDICT
pred_ntroglo <- pred_dist  
saveRDS(pred_ntroglo, "pred_ntroglo.rds")

#SALVANDO DADOS
dados_nt=dados_riqueza%>%
  dplyr::select(Riqueza,Distance)
dados_nt$Grupo <- "N_Troglo"
saveRDS(dados_nt, "dados_nt.rds")


