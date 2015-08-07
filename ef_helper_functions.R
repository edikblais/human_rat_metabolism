
# Factors can be very annoying so I typically set stringsAsFactors to FALSE
options(stringsAsFactors = FALSE)

# ef_df is a simple function used to clean up character and numeric columns in a data.frame
ef_df = function(x.tbl,x.filter = c(),
                 x.names = names(x.tbl)[sapply(x.tbl,function(x) {
                   !any(class(x) %in% c("numeric","integer","logical"))})]) {
  x.names = x.names[x.names %in% names(x.tbl)]
  if (length(x.names) < 1) return(x.tbl %>% as.tbl)
  x.col = x.names[1]
  x.tbl[[x.col]] = ef_df_cln(x.tbl[[x.col]])
  if (x.col %in% x.filter) {
    x.tbl = x.tbl[nchar(x.tbl[[x.col]])>0,]
  }
  ef_df(x.tbl,x.filter,x.names[-1]) %>% as.tbl
}

# ef_df_clean is called from ef_df to trim excess white space from characters
ef_df_cln = function(x) {
  if (any(class(x) %in% c("numeric","integer","logical"))) return(x)
  x = ifelse(!is.na(x),as.character(x),"")
  if (!any(grepl("\\s",x))) return(x)
  gsub("^\\s+|\\s+$","",x)
}

ef_reproduce_figure2 = function(gpr.info,gpr.limit = 9,gpr.legend = F) {
  gpr.color = c(`non-enzymatic` = "#616062",shared = "#674EA7",`rat-specific` = "#CC0000",`human-specific` = "#3C78D8")
  gpr.size = gpr.info %>% 
    mutate(rno = pmin(n_rno,gpr.limit),hsa = pmin(n_hsa,gpr.limit)) %>%
    mutate(organism = ifelse(rno > 0 & hsa > 0,"shared",ifelse(
      rno > 0,"rat-specific",ifelse(hsa > 0,"human-specific","non-enzymatic"))))
  gpr.size.count = gpr.size %>% count(rno,hsa) %>% ungroup %>% mutate(n_gpr = n) %>%
    mutate(organism = ifelse(rno > 0 & hsa > 0,"shared",ifelse(
      rno > 0,"rat-specific",ifelse(hsa > 0,"human-specific","non-enzymatic")))) %>% 
    mutate(color = gpr.color[organism])
  gpr.size.count %>% 
    ggplot(aes(x = factor(hsa), y = factor(rno)))+
    geom_tile(aes(fill = organism, alpha = log10(n_gpr)), width = 0.9, height = 0.9) + 
    scale_fill_manual(values = gpr.color) + scale_color_manual(values = gpr.color) +
    scale_alpha_continuous(range = c(.2,.9)) + 
    scale_x_discrete(breaks = 0:gpr.limit, labels = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit)) + 
    scale_y_discrete(breaks = 0:gpr.limit, labels = gsub(as.character(gpr.limit), paste0(gpr.limit,"+"), 0:gpr.limit)) + 
    xlab("Human GPR size") + ylab("Rat GPR size") + 
    theme_bw(base_size = 12) +
    theme(title = element_text(size = 12), legend.position = ifelse(gpr.legend, "top", "none"),
          axis.text.y = element_text(hjust = .7),
          panel.grid.major = element_line(size = NA),
          panel.grid.minor = element_line(size = NA))
}
ef_reproduce_figure3 = function(gpr.info,gpr.limit = 9) {
  gpr.color = c(`non-enzymatic` = "#616062",shared = "#674EA7",`rat-specific` = "#CC0000",`human-specific` = "#3C78D8")
  gpr.background = matrix(1,nrow = gpr.limit+1,ncol = gpr.limit+1,dimnames = list(
    0:gpr.limit,0:gpr.limit)) %>% melt %>% as.tbl %>% select(hsa = Var1,rno = Var2) %>% 
    mutate(organism = ifelse(rno > 0 & hsa > 0,"shared",ifelse(
      rno > 0,"rat-specific",ifelse(hsa > 0,"human-specific","non-enzymatic"))))
  gpr.size = gpr.info %>% mutate(rno = pmin(n_rno,gpr.limit),hsa = pmin(n_hsa,gpr.limit)) %>%
    mutate(organism = ifelse(rno > 0 & hsa > 0,"shared",ifelse(
      rno > 0,"rat-specific",ifelse(hsa > 0,"human-specific","non-enzymatic"))))
  gpr.size %>% ggplot(aes(x = factor(hsa), y = factor(rno)))+
    geom_point(aes(color = organism),alpha = 0.3, size = 1,
               position = position_jitter(width = 0.38, height = 0.38)) + #
    geom_tile(data = gpr.background,aes(fill = organism),
              color = NA,alpha = 0.1,width = .95,height = .95) +
    scale_fill_manual(values = gpr.color)+ scale_color_manual(values = gpr.color) +
    scale_x_discrete(breaks = 0:gpr.limit,labels = gsub(as.character(gpr.limit),paste0(gpr.limit,"+"),0:gpr.limit)) + 
    scale_y_discrete(breaks = 0:gpr.limit,labels = gsub(as.character(gpr.limit),paste0(gpr.limit,"+"),0:gpr.limit)) + 
    xlab("Human GPR size") + ylab("Rat GPR size") + 
    theme_bw(base_size = 12) +
    theme(legend.position="none",
          axis.text.y = element_text(hjust = .7),
          panel.grid.major = element_line(size = NA),
          panel.grid.minor = element_line(size = NA))
}






