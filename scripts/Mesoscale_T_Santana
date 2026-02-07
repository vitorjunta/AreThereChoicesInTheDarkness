#############################SANTANA GERAL###############################
#############################MESOESCALA################################

##################PREPARANDO O AMBIENTE R##############################
setwd("C:/Users/victor.junta/Desktop/SANTANA/ANALISES/MESOESCALA")
#setwd("G:/Meu Drive/BACKUP/UFLA/CEBS/Mestrado/Dissertação/Correção/GERAL/ANALISES/MESOESCALA")

library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(vegan)

##################TROGLOBITOS MICROESCALA##########################
##CARREGANDO OS DADOS##

#TABELA BASE
brutostroglo = read.csv("Geral_troglo_setores_r_limpo.csv", sep=";", header = T) %>%
  replace(is.na(.), 0)

brutostroglo$ID = paste(brutostroglo$Cave, brutostroglo$Sector, sep = "_")

#TABELA DE COMUNIDADE
comunidade_t = brutostroglo %>%
  select(ID, Clivina_sp1 : Xangoniscus._sp1)   

row.names(comunidade_t) = comunidade_t$ID
comunidade_t = comunidade_t[,-1]
write.table(comunidade_t, "Comunidade_T_Mesoescala.csv", row.names = T, sep = ";")

#TABELA DE COMUNIDADE BINARIA
comm.bin_t = decostand(comunidade_t, method = "pa")

#TABELA DE DADOS AMBIENTAIS
ambiente_t = brutostroglo %>%
  select(ID, Sector, Cave, Umidity, Temperature, Distance, Subs.Div., 
         Shelter.Div., Trophic.Div., Shelter.Av., Trophic.Av.)

row.names(ambiente_t) = ambiente_t$ID
ambiente_t = ambiente_t[,-1]

ambiente_t <- ambiente_t %>%
  mutate(across(where(is.numeric), scale))

#FILTRAR VARIAVEIS AMBIENTAIS
#tiff("ColinSGeral.tiff", units="in", width=20, height=20, res=300)
#chart.Correlation(ambiente_t[,-c(1,2, 3, 4)], histogram = TRUE, method = c("spearman"))
#dev.off()  
#remover subs.div.
#######DBRDA COM DISTANCIA DE JACCARD#### 

#FILTRANDO QUADRANTES ZERADOS
linhas_validas_t = rowSums(comm.bin_t) > 0

comm.bin_t = comm.bin_t[linhas_validas_t, ]
ambiente_t = ambiente_t[linhas_validas_t, ]

#MATRIX DE SIMILARIDADE DE JACCARD
diss_jaccard_t = vegdist(comm.bin_t, method = "jaccard")

#DBRDA
db0_t = capscale(diss_jaccard_t ~ 1, data = ambiente_t[,-c(1)])
db_t  = capscale(diss_jaccard_t ~ ., data = ambiente_t[,-c(1)])

summary(db_t)

#TESTES
set.seed(123)
anova(db_t)

set.seed(123)
anova(db_t, by = "term")

#ORDISTEP
set.seed(123)
modelo_final_t <- ordistep(db0_t, scope = formula(db_t), direction = "both")
formula(modelo_final_t)


#MODELO FINAL DBRDA 
db_final_t = capscale(diss_jaccard_t ~ Cave + Distance + Umidity, 
                      data = ambiente_t[,-c(1)])

summary(db_final_t)

set.seed(123)
anova(db_final_t)

set.seed(123)
anova(db_final_t, by="term")

#######GRAFICO DBRDA#### 
sites_scores_t = as.data.frame(db_final_t$CCA$u[, 1:2])
colnames(sites_scores_t) = c("CAP1", "CAP2")

env_vectors_t = as.data.frame(db_final_t$CCA$biplot[, 1:2])
colnames(env_vectors_t) = c("CAP1", "CAP2")
env_vectors_t$Variable = rownames(env_vectors_t)

sites_scores_t$Cave = ambiente_t$Cave
sites_scores_t$ID = rownames(sites_scores_t)

# Escalando vetores
env_vectors_t$Length = sqrt(env_vectors_t$CAP1^2 + env_vectors_t$CAP2^2)
scale_factor_t = 1.2 * max(abs(range(sites_scores_t[, 1:2]))) / max(env_vectors_t$Length)
env_vectors_t$CAP1_scaled = env_vectors_t$CAP1 * scale_factor_t
env_vectors_t$CAP2_scaled = env_vectors_t$CAP2 * scale_factor_t

# Variância explicada
var_cap1_t = round(100 * db_final_t$CCA$eig[1]/sum(db_t$CCA$eig), 1)
var_cap2_t = round(100 * db_final_t$CCA$eig[2]/sum(db_t$CCA$eig), 1)
var_total_t = round(100 * db_final_t$CCA$tot.chi/db_t$tot.chi, 1)

#nomes nonitos
sites_scores_t <- sites_scores_t %>%
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
  "Cedrão Cave"                = "#ff7f0e",   # Verde (mantida da original)
  "Cedro Cave"                 = "#ff7f0e",   # Verde (mantida da original)
  "Pedra Escrevida Cave"       = "#9467bd",   # Vermelho
  "Labirinto Cave"             = "#9467bd",   # Roxo
  "Leão Cave"                  = "#ff7f0e",   # Marrom (mantida da original)
  "Padre Cave"                 = "#ff7f0e",   # Laranja (mantida da original)
  "Racha Bovina Cave"          = "#9467bd",   # Azul (mantida da original)
  "Boqueirão Cave"             = "#9467bd",   # Verde-amarelado
  "Couve Flor Cave"            = "#9467bd",   # Azul-petróleo
  "São Geraldo Cave"           = "#9467bd",   # Laranja-claro
  "Tunel 2 Cave"               = "#9467bd",   # Verde-claro
  "Geraldo Cruz Cave"          = "#9467bd",   # Vermelho-claro
  "Duas Cobras Cave"           = "#9467bd",   # Lilás
  "Olho D'Água do Cumbra Cave" = "#9467bd",   # Bege
  "Fenda Oblíquoa Cave"        = "#9467bd",   # Rosa-claro
  "Tunel 1 Cave"               = "#9467bd",   # Amarelo-claro
  "Pajeú Cave"                 = "#9467bd",   # Azul-claro
  "Grota Cave"                 = "#ff7f0e",   # Azul-escuro
  "Salobro Cave"               = "#9467bd",   # Verde-escuro
  "Cedrículo Cave"             = "#9467bd",   # Marrom-claro
  "Cinquentona Cave"           = "#9467bd",   # Vinho
  "Escrevidinha Cave"          = "#ff7f0e",   # Roxo-escuro
  "Cristal Cave"               = "#9467bd"    # Azul-claro pastel
)


# DEFINIR POSIÇÕES MANUAIS PARA OS NOMES DAS CAVERNAS
posicoes_nomes <- data.frame(
  Cave_formatted = c(
    "Boqueirão Cave",
    "Cânion da Baixa Verde Cave",
    "Cedrão Cave",
    "Cedro Cave",
    "Escrevidinha Cave",
    "Labirinto Cave",
    "Leão Cave",
    "Padre Cave",
    "Pedra Escrevida Cave",
    "Racha Bovina Cave",
    "Salobro Cave",
    "São Geraldo Cave",
    "Tunel 2 Cave"
  ),
  x_pos = c(0.106, 0.198, -0.0597, -0.135,  -0.108,
            0.132,  -0.0865,  -0.0817, 0.169, 0.218,
            0.0376,  0.208, 0.208),
  y_pos = c(-0.00429, 0.00636, -0.0388, -0.0950,  0.267,
            0.00477, -0.00465,  -0.0388, 0.000653, 0.00636,
            0.02263,   0.00636,   0.02638)
) %>%
  dplyr::mutate(ID = case_when(
    Cave_formatted == "Cânion da Baixa Verde Cave" ~ 1,
    Cave_formatted == "Padre Cave"                 ~ 2,
    Cave_formatted == "Labirinto Cave"             ~ 3,
    Cave_formatted == "Boqueirão Cave"             ~ 4,
    Cave_formatted == "Pedra Escrevida Cave"       ~ 5,
    Cave_formatted == "Tunel 2 Cave"               ~ 7,
    Cave_formatted == "Racha Bovina Cave"          ~ 10,
    Cave_formatted == "Cedro Cave"                 ~ 15,
    Cave_formatted == "Cedrão Cave"                ~ 16,
    Cave_formatted == "Salobro Cave"               ~ 20,
    Cave_formatted == "Leão Cave"                  ~ 22,
    Cave_formatted == "Cinquentona Cave"           ~ 23,
    Cave_formatted == "Escrevidinha Cave"          ~ 24,
    TRUE ~ NA_integer_
  ))


# Gráfico
ggplot() +
  #stat_ellipse(data = sites_scores_t, 
   #            aes(x = CAP1, y = CAP2, color = Cave_formatted),
    #           type = "norm", level = 0.95, linewidth = 0.8, alpha = 0.6) +
  geom_point(data = sites_scores_t, 
             aes(x = CAP1, y = CAP2, color = Cave_formatted),
             size = 6, alpha = 0.7) +
  geom_segment(data = env_vectors_t[env_vectors_t$Variable == "Distance",],
               aes(x = 0, y = 0, xend = CAP1_scaled, yend = CAP2_scaled),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "black", linewidth = 1.2) +
  geom_text_repel(data = env_vectors_t[env_vectors_t$Variable == "Distance",],
                  aes(x = CAP1_scaled * 1.1, y = CAP2_scaled * 1.1, 
                      label = Variable),
                  color = "black", fontface = "bold", size = 6) +
  geom_text(data = posicoes_nomes,
            aes(x = x_pos, y = y_pos, label = ID, color = Cave_formatted),
            size = 6, alpha = 1, fontface = "bold", check_overlap = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
  labs(title = "Mesoscale - dbRDA - Troglobitic species",
       subtitle = paste("Jaccard distance | ", var_total_t, "% explained"),
       x = paste("CAP1 (", var_cap1_t, "%)"),
       y = paste("CAP2 (", var_cap2_t, "%)")) +
  scale_color_manual(values = cores_cavernas) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.9),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(face = "bold", size = 16),
        legend.position = "none")


ggsave("DBRDA_Troglos_Meso.tiff", width = 14, height = 10, dpi = 300, bg = "white")

###################################################################################
#########################ENTENDENDO AS RELAÇÕES ESPACIAIS#####################
brutostroglo$ID


coord_setor_t = read.csv("setores_utm_ID.csv", sep=";", header = T)
coord_setor_t <- coord_setor_t[coord_setor_t$ID %in% rownames(ambiente_t), ]
coord_setor_t$ID

coords_matrix_t = cbind(coord_setor_t$x, coord_setor_t$y)
dist_geo_t = dist(coords_matrix_t)

mantel_result_t = mantel(diss_jaccard_t, dist_geo_t, 
                          method = "spearman",
                          permutations = 9999)

distancias_t = as.vector(as.matrix(dist_geo_t))
dissimilaridades_t = as.vector(as.matrix(diss_jaccard_t))

dados_plot_t = data.frame(
  Distancia_m = distancias_t,
  Distancia_km = distancias_t/1000,
  Dissimilaridade = dissimilaridades_t
)

ggplot(dados_plot_t, aes(x = Distancia_km, y = Dissimilaridade)) +
  geom_point(alpha = 0.2, size = 2, color = "gray40") +
  geom_smooth(method = "gam", color = "blue", se = TRUE, linewidth = 1.2) +
  labs(
    title = "Geographic Distance vs. Dissimilarity",
    subtitle = paste("Mantel r =", round(mantel_result_t$statistic, 4),
                     "| p =", mantel_result_t$signif),
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


ggsave("Mantel_Troglos_Setor.tiff", width = 14, height = 10, dpi = 300, bg = "white")


######################PARTIÇÃO DA VARIAÇÃO###################################
###DADOS AMBIENAIS###
X_env = ambiente_t[, c("Distance")]

###GERANDO POLINOMIOS DE TENDENCIA ESPACIAL###

coords = coord_setor_t[, c("x", "y")]
rownames(coords) = coord_setor_t$ID

dist_geo = dist(coords)

library(adespatial)
dbmem = dbmem(dist_geo, silent = FALSE) 

X_spa = as.data.frame(dbmem)

#PARTIÇÃO DA VARIANCIA
vp = varpart(comm.bin_t, ~ Cave + Distance, # Environment
             X_spa, 
             data = ambiente_t)

vp
plot(vp,
     digits = 2,
     bg = c("skyblue", "pink"),
     Xnames = c("Ambiente", "Espaço"))

#SIGNIFICANCIA DA FRAÇÃO PURAMENTE AMBIENTAL
rda_env = rda(comm.bin_t ~ Cave + Distance + Condition(MEM1 + MEM2),
              data = cbind(ambiente_t, X_spa))

set.seed(125)
anova(rda_env)

set.seed(125)
anova(rda_env, by = "margin")

#SIGNIFICANCIA DA FRAÇÃO PURAMENTE ESPACIAL
rda_spa = rda(comm.bin_t ~ MEM1 + MEM2  + Condition(Cave + Distance),
              data = cbind(ambiente_t, X_spa))

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

r2_shared = max(0, r2_env + r2_spa - r2_total)
exclusive_env = max(0, r2_env - r2_shared)
exclusive_spa = max(0, r2_spa - r2_shared)

fit = euler(c(
  "Environmental" = exclusive_env,
  "Spatial" = exclusive_spa,
  "Environmental&Spatial" = r2_shared
))


tiff("Venn_Meso_T.tiff", 
     width = 16, height = 12, units = "cm", res = 600, compression = "lzw")

plot(fit,
     quantities = list(type = "percent", fontsize = 10),
     fills = list(fill = c("#1f77b4", "#ff7f0e", "#2ca02c"), alpha = 0.7),
     labels = list(font = 2),
     main = "Variance Partitioning - Troglos")

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
brutostroglo = read.csv("Geral_troglo_setores_r_limpo.csv", sep=";", header = T)%>%
  replace(is.na(.), 0)

brutostroglo$ID = paste(brutostroglo$Cave, brutostroglo$Sector, sep = "_")

#TABELA DE COMUNIDADE
comunidade=brutostroglo%>%
  dplyr::select(ID, Clivina_sp1 : Xangoniscus._sp1)

row.names(comunidade) = comunidade$ID
comunidade=comunidade[,-1]

#TABELA DE DADOS AMBIENTAIS
ambiente=brutostroglo%>%
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
pred_dist = ggpredict(m1, terms = "Umidity")
plot(pred_dist) + labs(x = "Humidity", y = "Predicted richness")

#SALVANDO PREDICT
pred_troglo <- pred_dist  
saveRDS(pred_troglo, "pred_troglo.rds")

#SALVANDO DADOS
dados_t=dados_riqueza%>%
  dplyr::select(Riqueza,Distance)
dados_t$Grupo <- "Troglo"
saveRDS(dados_t, "dados_t.rds")

