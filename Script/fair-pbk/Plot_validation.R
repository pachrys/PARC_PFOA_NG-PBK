library(readr)

load(here(OUTPUT, "RawData_PFOA_PBK.RData"))
R <- select(RawData$PBK_OUTPUT, c(time, !starts_with("C"))) %>% 
  select(!c(OD, Ain, Atot, MB))
R <- R %>% pivot_longer(AIL:AMp, names_to = "Tissue", values_to = "Amount") %>% 
  mutate(language= "R")

SBML <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Output/oral/oral_single_NG_PBK_PFOA.csv")

replacements <- c(
  "AIntestine_lumen" = "AIL",
  "AIntestine_tissue" = "AI",
  "AFeces" = "AFe",
  "ALiver_extracellular" = "AL_ec",
  "ALiver_intracellular" = "AL_ic",
  "AKidney_proximal_tubule_tissue" = "APTT",
  "AKidney_proximal_tubule_lumen" = "APTL",
  "AKidney_rest_tissue" = "ARKT",
  "AKidney_rest_lumen" = "ARKL",
  "AUrine" = "AUr",
  "AAdipose" = "AA",
  "ARest" = "AR",
  "APlasma" = "AP",
  "AMenstrual" = "AMp"
)

SBML <- SBML %>% set_names(~(.) %>%
            str_replace_all(replacements))

SBML <- SBML %>% pivot_longer(AIL:AMp, names_to = "Tissue", values_to = "Amount") %>% 
  mutate(language= "SBML")

results <- rbind(R, SBML)
write.csv(results, here(OUTPUT,
                        "oral_single_R_SBML_results.csv"),
          row.names = FALSE)

plot <- results %>% 
  ggplot(aes(time, Amount, colour = language)) +
  geom_path() +
  facet_wrap(~Tissue, scales = "free") +
  theme_minimal()
plot
ggsave(
  filename = here(OUTPUT,"oral_single_R_SBML_result.png")
)