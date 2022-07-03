angle <- function(x,y){
  dot.prod <- x %*% y 
  norm.x <- norm(x, type = "2")
  norm.y <- norm(y, type = "2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}

full.cube <- data.frame(data.table::fread('F_me.csv', fill = T))
full.cube <- data.frame(apply(full.cube, 2, as.numeric))
full.cube <- full.cube[-(1:2),]
n.atoms <- as.numeric(full.cube[1,1])
structure <- full.cube[1:(4 + n.atoms), 1:5]
density <- full.cube[-c(1:(4 + n.atoms)),]
x.origin <- structure[1, 2]
y.origin <- structure[1, 3]
z.origin <- structure[1, 4]
x.steps <- structure[2, 1]
y.steps <- structure[3, 1]
z.steps <- structure[4, 1]
x.size <- structure[2, 2]
y.size <- structure[3, 3]
z.size <- structure[4, 4]
xyz <- structure[5:(4 + n.atoms), c(1, 3:5)]
names(xyz) <- c('atom', 'x', 'y', 'z')
suppressMessages(xyz$atom <- plyr::mapvalues(xyz$atom,
                                             from = c("1", "5", "6", "7", "8", "9", "14", "15", "16", "17", "35", "53","27",'28'),
                                             to = c(
                                               "H", "B", "C", "N", "O", "F", "Si", "P",
                                               "S", "Cl", "Br", "I","Co", "Ni"
                                             )
))

blocks <- seq(1, nrow(density), ceiling(z.steps/6))
row.dens <- data.frame(matrix(nrow = length(blocks),
                              ncol = ceiling(pracma::nthroot(x.steps * y.steps * z.steps, 3))))

for (i in blocks) {
  vec <- c(t(density[i:(i + (ceiling(z.steps/6) - 1)), ]))
  vec[which(is.na(vec))] <- 0
  row.dens[which(blocks == i), ] <- vec
}
for (i in 1:nrow(row.dens)) {
  for (j in 1:ncol(row.dens)) {
    if (row.dens[i, j] < 0.002) row.dens[i, j] <- 0
  }
}

x.blocks <- list()
for (i in seq(1, nrow(row.dens), y.steps)) {
  x.blocks[[i]] <- row.dens[i:(i + (y.steps - 1)), ]
}
x.blocks <- x.blocks[seq(1, nrow(row.dens), y.steps)]
point.space <- x.blocks

for (block in 1:length(x.blocks)) {
  for (row in 1:nrow(x.blocks[[block]])) { 
    for (col in 1:ncol(x.blocks[[block]][row, ])) {
      if (point.space[[block]][row, col] != 0) {
        point.space[[block]][row, col] <- paste(as.character(block - 1),
                                                as.character(row - 1),
                                                as.character(col - 1)) 
      }
    }
  }
}

dense.points <- vector()
for (i in 1:length(point.space)) {
  dense.points <- append(dense.points, unlist(point.space[[i]])[which(point.space[[i]] != 0)])
}
coordinate.space <- matrix(nrow = length(dense.points), ncol = 3)
colnames(coordinate.space) <- c('x','y','z')

for (i in 1:nrow(coordinate.space)) {
  num.pt <- as.numeric(unlist(strsplit(dense.points[i], " ")))
  coordinate.space[i, 1] <- x.origin + num.pt[1] * x.size
  coordinate.space[i, 2] <- y.origin + num.pt[2] * y.size
  coordinate.space[i, 3] <- z.origin + num.pt[3] * z.size
}
rgl::plot3d(coordinate.space)
coordinate.space <- rbind(xyz[,2:4], coordinate.space)

atoms <- strsplit(coor.atoms, " ")
unlisted.atoms <- unlist(atoms)
numeric.atoms <- as.numeric(unlisted.atoms)
xyz <- coordinate.space
mag <- function(vector) {
  sqrt(vector[[1]]^2 + vector[[2]]^2 + vector[[3]]^2)
}
new_origin <- xyz[numeric.atoms[[1]], 1:3]
new_y <- as.numeric((xyz[numeric.atoms[[2]], 1:3] - new_origin) /
                      mag(xyz[numeric.atoms[[2]], 1:3] - new_origin))
coplane <- as.numeric((xyz[numeric.atoms[[3]], 1:3] - new_origin) /
                        mag(xyz[numeric.atoms[[3]], 1:3] - new_origin))
cross_y_coplane <- pracma::cross(coplane, new_y)
coef_mat <- aperm(array(c(
  new_y,
  coplane,
  cross_y_coplane
),
dim = c(3, 3)
))
angle_new.y_coplane <- angle(coplane, new_y)
x_ang_new.y <- pi / 2
cop_ang_x <- angle_new.y_coplane - x_ang_new.y
result_vec <- c(0, cos(cop_ang_x), 0)
new_x <- solve(coef_mat, result_vec)
new_z <- pracma::cross(new_x, new_y)
new_basis <- aperm(array(c(new_x, new_y, new_z), dim = c(3, 3)))
new_origin <- xyz[numeric.atoms[[1]], 1:3]
new_coordinates <- matrix(nrow = dim(xyz)[[1]], ncol = 3)
transformed_coordinates <- matrix(nrow = dim(xyz)[[1]], ncol = 3)
for (i in 1:dim(xyz)[[1]]) {
  new_coordinates[i, ] <- as.numeric(xyz[i, 1:3] - new_origin)
  transformed_coordinates[i, ] <- aperm(new_basis %*% new_coordinates[i, ])
}
transformed_coordinates <- round(transformed_coordinates, 6)
rgl::plot3d(transformed_coordinates)


  
  
