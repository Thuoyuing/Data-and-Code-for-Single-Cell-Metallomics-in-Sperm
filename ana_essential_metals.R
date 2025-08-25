## library packages ---------------------------------------------------------
if(!require(pacman)){
  install.packages("pacman")
}
pacman::p_load(tidyverse,
               magrittr,
               mlr3verse,
               DALEX,
               data.table,
               lme4,
               lmerTest,
               caret,
               Rcpp)
rm(list = ls())

## define data path -----------------------------------------------------------
cell_data_path <- "data/cell_data/"
outcome_data <- "data/outcome_data/outcdt.csv"

## ReadMT -----------------------------------------------------------
source("functions/ReadMT.R")
ReadMT(cell_data_path,
       outcome_data) -> MT

MT$Dir.RawData$dir.outc
MT$Dir.RawData$dir.cell

MT$RawData$dt_outc
MT$RawData$dt_cell


## load cell data -----------------------------------------------------------
# Not recommended
# Especially when there is insufficient memory to read large-scale single-cell data.
source("functions/LoadSCDT.R")
LoadSCDT(MT) -> MT

MT$RawData$dt_cell %>% names()

## data reshape -----------------------------------------------------------
source("functions/DtReshape.R")
source("functions/PerctCal.R")
source("functions/QuanCal.R")
sourceCpp("functions/HCalC.cpp")


DtReshape(MT,
          num_quan = 100,
          q.cut = 0,
          heterog = T,
          prevelen = T) -> MT


MT$ReshapedDT
MT$RawData$dt_outc

MT$ReshapedDT %>% names() %>% 
  tibble(value=.) %>% 
  dplyr::filter(stringr::str_detect(value,"_P")) %>% 
  as.matrix() %>% as.vector() -> vars.p
MT$ReshapedDT %>% 
  dplyr::select(all_of(vars.p)) %>% 
  summarise(across(where(is.numeric), max, na.rm = TRUE)) %>% 
  t()
  
# It is recommended to save the RData at this step, after prolonged processing
save(MT,file = "MT.RData")

load("MT.RData")
## data bind -----------------------------------------------------------
source("functions/DtBind.R")
DtBind(MT) -> MT

MT$ReshapedDT %>% names()

## data summary -----------------------------------------------------------
source("functions/FtSmry.R")
FtSmry(MT,SC=T) -> MT
MT$FeaturePerformance$plot
MT$FeaturePerformance$table %>% 
  ggplot(aes(x=FeatureType,y=remain,fill=FeatureType))+
  geom_col(width = 0.75)+
  scale_fill_manual(values=c("Q"="#D75615",
                             "H"="#EDB11A",
                             "M"="#FF6B6B",
                             "P"="#7E3E8A",
                             "SC"="#78AB31"))+
  labs(x="Feature type",
       y="Non-zero-variance \nelement proportion")+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))+
  mytheme
ggsave("output/plot/NonZeroVar.png",
       width = 4,
       height = 4)
ggsave("output/plot/NonZeroVar.svg",
       width = 4,
       height = 4)
## feature trans -----------------------------------------------------------
source("functions/FtTrans.R")
FtTrans(MT,
        features=c("Q","P","H","Mean"),
        variable=c("motility"),
        logtrans=T,
        scale=T) -> MT

MT$TransReshapedDT
save(MT,file = "MT.RData")

load("MT.RData")
## feature evaluation -----------------------------------------------------------
source("functions/FtEvalu.R")
FtEvalu(MT,
        features=c("Q","P","H","Mean"),
        outcome="motility",
        covariate=c("age"),
        SC=T,
        task="20250529_TXW") -> MT

MT$RawData$dt_cell %>% 
  mutate(filename=stringr::str_extract(dir.path,"(?<=/)[^/]+(?=\\.)")) %>% 
  left_join(.,
            MT$RawData$dt_outc,
            by="filename") -> dt.sc
dt.sc %>% 
  dplyr::select(-all_of(c("filename","dir.path"))) %>% 
  names()

MT$RawData$dt_cell %>% 
  dplyr::select(-dir.path) %>% 
  names() -> Xs

Cova <- "age"

dt.sc %>% 
  mutate(across(all_of(Xs),  
                ~log(.+0.001))) %>% 
  mutate(across(all_of(Xs),  
                ~scale(.))) -> dt.sc


dt.sc %>%
  dplyr::select(Xs) %>%
  caret::nearZeroVar(saveMetrics = T,
                     names = T) %>%
  dplyr::filter(nzv == "TRUE") %>%
  rownames()  -> varsdeleted_sc

Xs <- setdiff(Xs,varsdeleted_sc)

tibble(Y="motility",
       X=Xs) %>% 
  mutate(Formula=stringr::str_c(Y,"~",X,
                                "+(1|filename)+",Cova)) %>% 
  mutate( Model = map(Formula,
                      ~lmer(as.formula(.),
                            data=dt.sc)),
          Coef.dt = map(Model,
                        ~as.data.frame(summary(.)$coefficients) %>% set_names(str_c("",names(.)))),
          Coef.dt = map2(Coef.dt,
                         Model,
                         ~inner_join(bind_cols(Variable = rownames(.x),.x),
                                     bind_cols(Variable = rownames(confint.merMod(.y,method = "Wald")),
                                               low = confint.merMod(.y,method = "Wald")[,"2.5 %"],
                                               high = confint.merMod(.y,method = "Wald")[,"97.5 %"]),
                                     by = "Variable"))) %>%
  unnest(Coef.dt) %>%
  mutate(N_Obs = map_dbl(Model,
                         ~length(summary(.)$residuals)),
         N_Subject = map_dbl(Model,
                             ~summary(.)$ngrps)) %>% 
  dplyr::select(-Model) %>% 
  dplyr::filter(X==Variable) -> res.sc
res.sc
res.sc %>% 
  dplyr::filter(`Pr(>|t|)`<0.05) %>% 
  nrow() -> sig.sc

MT$TransReshapedDT %>% 
  names() %>% tibble(value=.) %>% 
  dplyr::filter(stringr::str_detect(value,"(?<=_)[A-Z]")) %>% 
  as.matrix() %>% as.vector() -> vars.features

MT$TransReshapedDT %>%
  dplyr::select(vars.features) %>%
  caret::nearZeroVar(saveMetrics = T,
                     names = T) %>%
  dplyr::filter(nzv == "TRUE") %>%
  rownames()  -> varsdeleted_features

vars.features <- setdiff(vars.features,varsdeleted_features)


tibble(Y="motility",
       X=vars.features) %>% 
  mutate(Formula=map2_chr(X,Y,
                          ~stringr::str_c(.y,"~",.x,"+",Cova))) %>% 
  mutate(Model = map(Formula,
                     ~lm(as.formula(.),
                         data = MT$TransReshapedDT))) %>% 
  mutate(Coef = map(Model,
                    ~summary(.) %>% 
                      coef() %>% 
                      as.data.frame %>% 
                      rownames_to_column() %>% 
                      as_tibble() %>% rename("Vars" = rowname)),
         Confint = map(Model,
                       ~confint.lm(.,level = 0.95) %>% 
                         as.data.frame() %>% 
                         rownames_to_column() %>% set_names(c("Vars","Low","High"))),
         table = map2(Coef,
                      Confint,
                      ~inner_join(.x,
                                  .y,
                                  by = "Vars"))) %>% 
  unnest(table) %>% 
  mutate(EXP_Value = map_dbl(Estimate,exp),
         EXP_Low = map_dbl(Low,exp),
         EXP_High = map_dbl(High,exp),
         N = map_dbl(Model,
                     ~(length(.$fitted.values)))) %>% 
  dplyr::filter(Vars != "(Intercept)") %>% 
  dplyr::select(Y,
                X,
                Formula,
                Vars,
                Estimate,
                Low,
                High,
                "Pr(>|t|)",
                EXP_Value,
                EXP_Low,
                EXP_High) %>% 
  set_names(c("Y","X","Formula","Vars","Beta","Low",
              "High","P.value","EXP.Beta","EXP.Low","EXP.High")) %>% 
  dplyr::filter(X==Vars) -> res.features

res.features %>% 
  dplyr::filter(P.value<0.05) %>% 
  dplyr::select(X) %>% 
  mutate(Feature=stringr::str_extract(X,"(?<=_).+")) %>% 
  mutate(Feature=ifelse(stringr::str_detect(Feature,"Q[0-9].+"),"Quantile",Feature)) %>% 
  mutate(Element=stringr::str_extract(X,".+(?=_)")) -> lm.sig.tab


theme(text = element_text(family = 'sans'),     # 字体
      plot.title = element_text(size='17',
                                hjust = 0.5),   
      axis.title= element_text(size='14'),     
      axis.text= element_text(size='12'),
      legend.title = element_text(size=14),
      legend.text = element_text(size=12),
      strip.text = element_text(size = 14)) -> mytheme 
lm.sig.tab %>% 
  dplyr::select(Feature) %>% 
  table %>% t %>% 
  as.data.frame() %>% 
  dplyr::select(-Var1) %>% 
  add_row(Feature="SC",
          Freq=sig.sc) %>% 
  ggplot(aes(x=Feature,y=Freq,colour = Feature,fill = Feature))+
  geom_col(width = 0.75)+
  theme_classic()+
  # scale_y_break(c(30, 220),
  #               space = 0.25) +
  scale_colour_manual(values=c("Quantile"="#D75615",
                               "H"="#EDB11A",
                               "Mean"="#FF6B6B",
                               "P"="#7E3E8A",
                               "SC"="#78AB31",
                               "H"="#78AB31"))+
  scale_fill_manual(values=c("Quantile"="#D75615",
                             "H"="#EDB11A",
                             "Mean"="#FF6B6B",
                             "P"="#7E3E8A",
                             "SC"="#78AB31",
                             "H"="#78AB31"))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))+
  labs(x="Feature type",
       y="Number of \nsignificant featurs")+
  mytheme+
  theme(legend.position = "none")
ggsave(filename = stringr::str_c("output/plot/sigcol.png"),
       width = 4,
       height = 4)
ggsave(filename = stringr::str_c("output/plot/sigcol.svg"),
       width = 4,
       height = 4)
## feature identification -----------------------------------------------------------
source("functions/FtIdent.R")
FtEvalu(MT,
        features=c("Q","P","H","Mean"),
        outcome="motility",
        covariate=c("age"),
        task="20250529_TXW") -> MT
