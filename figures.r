library(gridExtra)
library(ggplot2)
library(tidyr)

Figure2Adata <- read.csv("Figure2Adata.csv")[,-1]
Figure2Bdata <- read.csv("Figure2Bdata.csv")[,-1]

# Reformat 2a
long2a <- gather(Figure2Adata, key="key", value="value", -MaxPleioSize)
long2a$label <- "Cochran's Q"
long2a$label[long2a$key == "PowerMREgg"] <- "MR-Egger intercept"
long2a$label[long2a$key == "PowerRuckQ"] <- "Rucker's Q"

# Reformat 2b
Figure2Bdata$SNP <- 1:25
thresh <- data.frame(thresh=qchisq(c(0.95, 1-0.05/25),1), lab=c("P = 0.05", "P = 0.05/25"))
long2b <- gather(Figure2Bdata, key="key", value="value", -SNP)
long2b$label <- "Cochran's Q"
long2b$label[long2b$key == "RuckerQ_contrib"] <- "Rucker's Q"

# Draw 2a
p1 <- ggplot(long2a, aes(x=MaxPleioSize, y=value, group=label)) +
	geom_line(aes(colour=label, group=label), size=1) +
	geom_hline(yintercept=0.05, linetype="dotted") +
	scale_color_brewer(type="qual") +
	labs(x="Maximum Size of Direct/Pleiotropic Effect", y="Power to Detect Global Pleiotropy", colour=NULL, title="A") +
	theme(
		legend.position=c(0.2, 0.85),
		legend.box.background = element_rect(size=0.2, fill=NULL),
		legend.background = element_rect(fill = alpha("red", 0)),
		axis.line=element_line(size=0.2),
		panel.grid.minor=element_blank(),
		panel.grid.major=element_blank(),
		panel.background = element_rect(fill = alpha("red", 0)),
		plot.title=element_text(hjust=-0.15, size=25),
		axis.title=element_text(size=10, family="Helvetica"),
		axis.title.x = element_text(margin=margin(t=20, r=0, b=0, l=0)),
		axis.title.y = element_text(margin=margin(t=0, r=20, b=0, l=0)),
		axis.text=element_text(size=10, family="Helvetica", colour="black"),
		legend.key=element_blank()
	)
p1
ggsave(p1, file="figure2a.pdf", width=6, height=6)
ggsave(p1, file="figure2a.eps", width=6, height=6)


# Draw 2b
p2 <- ggplot(long2b, aes(x=SNP, y=value)) +
	geom_hline(data=thresh, aes(yintercept=thresh, linetype=lab)) +
	geom_point(aes(shape=label), size=3) +
	labs(x="SNP", y="Q Statistic Contribution", linetype="Threshold", shape="Method", title="B") +
	# theme_bw() +
	scale_shape_manual(values=c(2,16)) +
	scale_linetype_manual(values=c("dotted", "dashed")) +
	theme(
		legend.position=c(0.2, 0.8),legend.box.background = element_rect(size=0.2, fill=NULL),legend.background = element_rect(fill = alpha("red", 0)), axis.line=element_line(size=0.2),
		panel.grid.minor=element_blank(), panel.grid.major=element_blank(),panel.background = element_rect(fill = alpha("red", 0)),
		plot.title=element_text(hjust=-0.15, size=25),
		axis.title=element_text(size=10, family="Helvetica"),
		axis.title.x = element_text(margin=margin(t=20, r=0, b=0, l=0)),
		axis.title.y = element_text(margin=margin(t=0, r=20, b=0, l=0)),
		axis.text=element_text(size=10, family="Helvetica", colour="black"), legend.key=element_blank(), legend.title=element_blank(),
		legend.spacing=unit(-0.5, "cm")
	) +	ylim(0, 20)
p2
ggsave(p2, file="figure2b.pdf", width=6, height=6)
ggsave(p2, file="figure2b.eps", width=6, height=6)

# Plot
pdf("figure2.pdf", width=12, height=6)
grid.arrange(p1, p2, ncol=2)
dev.off()

setEPS()
postscript("figure2.eps", width=12, height=6)
grid.arrange(p1, p2, ncol=2)
dev.off()
