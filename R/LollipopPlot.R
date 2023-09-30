#' Plotting ParsimonyC3 results with Lollipop Plot
#'
#' @param The results from ParsimonyC3
#' @param top The top number of LR combinations, default: 10
#'
#' @return Plotting the priortized LR pair combinations
#' @examples
#' \donttest{
#' LollipopPlot(result,top = 10)
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @export

LollipopPlot = function(
    df,
    top=10
) {
  df_s=df %>% top_n(n = -top, wt = z_shortest)
  df_s=df_s[order(-df_s$z_shortest),]
  df_s$num=paste0("C",1:top)
  df_s$pairs=factor(df_s$pairs,levels = df_s$pairs)
  p1=ggplot(df_s,aes(z_shortest,pairs,fill=num))+
    geom_segment(aes(yend=pairs,xend=0),size=0.8,color="black")+
    geom_point(aes(size=3),pch=21)+
    labs(x=NULL,y=NULL)+
    ggtitle(label="shortest.paths")+
    xlab("Z score")+
    theme(legend.position = "none",
          plot.title=element_text(size =14,color="black",vjust=0.5,hjust=0.5),
          axis.text.y = element_text(size =12,angle = 0, vjust = 0.5,color="black"),
          axis.text.x = element_text(size =12, angle = 0,color="black"),
          axis.title = element_text(size=12),
          plot.margin = unit(c(1,1,1,1),units="cm"),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x.top = element_blank(),
          axis.ticks.x.top = element_blank())+
    scale_x_continuous(guide = guide_axis(angle=0))
  df_f=df %>% top_n(n = -top, wt = z_diffusion)
  df_f=df_f[order(-df_f$z_diffusion),]
  df_f$num=paste0("C",1:top)
  df_f$pairs=factor(df_f$pairs,levels = df_f$pairs)
  p2=ggplot(df_f,aes(z_diffusion,pairs,fill=num))+
    geom_segment(aes(yend=pairs,xend=0),size=0.8,color="black")+
    geom_point(aes(size=3),pch=21)+
    labs(x=NULL,y=NULL)+
    ggtitle(label="diffusion")+
    xlab("Z score")+
    theme(legend.position = "none",
          plot.title=element_text(size =14,color="black",vjust=0.5,hjust=0.5),
          axis.text.y = element_text(size =12,angle = 0, vjust = 0.5,color="black"),
          axis.text.x = element_text(size =12, angle = 0,color="black"),
          axis.title = element_text(size=12),
          plot.margin = unit(c(1,1,1,1),units="cm"),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x.top = element_blank(),
          axis.ticks.x.top = element_blank())+
    scale_x_continuous(guide = guide_axis(angle=0))
  p1|p2
}
