mean_and_var <- function(Wline.mat,Zline.mat) {
x_1 <- colSums(Wline.mat * Zline.mat)
x_2 <- colSums(Wline.mat * Zline.mat^2)
v <- x_2 - (x_1)^2
list("means"=x_1,"vars"=v)
}