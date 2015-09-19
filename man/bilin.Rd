\name{bilin}
\alias{bilin}
\title{Fast 2D Linear interpolation}
\usage{
bilin( pt, x, y, f, n_x, n_y)
}
\arguments{
  \item{pt}{An \eqn{(x,y)} vector at which to approximately evaluate a function}
  \item{x}{A vector of x-values}
  \item{y}{A vector of y-values}
  \item{f}{A vector of z-values ordered first by \eqn{x} and then by \eqn{y}.  If \eqn{x=(x_1,x_2)} and \eqn{y=(y_1,y_2)} then \eqn{f = (f(x_1,y_1),f(x_1,y_2),f(x_2,y_1),f(x_2,y_2))}}
  \item{n_x}{An integer.  The length of \code{x}}
  \item{n_y}{An integer.  The length of \code{y}}
}
\value{
Returns the bilinear interpolation of \eqn{(x,y,f)} evaluated at \code{pt}}.
}
\description{
  Provides very fast 2D linear interpolation.
}
\examples{
## 2D Interpolation
f <- function(pt) pt[1] * cos(pt[2]) + pt[2] ^ 2 * sin(pt[1])
    # An arbitrary function
n.coarse <- 7
    # Number of points at which the interpolating grid is defined
n.fine <- 201
    # Number of poiints to evaluate the interpolation
xy.coarse <- cbind( rep( 1:n.coarse, each=n.coarse), rep( 1:n.coarse, n.coarse ) )
    # The matrix of coarse (x,y) values at which to evaluate f
v.f.w <- apply( xy.coarse, 1, f )
m.f.w <- matrix( v.f.w, n.coarse, n.coarse )
    # The vector & matrix form of the z-values

xy.fine <- seq( 1, n.coarse, length.out=n.fine )
    # A fine grid on which to evaluate the interpolation
m.xy.fine <- cbind( rep( xy.fine, each=n.fine), rep( xy.fine, n.fine ) )
    # The matrix form
v.f.v <- apply( m.xy.fine, 1, bilin, x=1:n.coarse, y=1:n.coarse, f=v.f.w, n_x=n.coarse, n_y=n.coarse )
    # The linear interpolation
v.f.true <- apply( m.xy.fine, 1, f )
    # The true value
m.f.v <- matrix( v.f.v, n.fine, n.fine )
m.f.true <- matrix( v.f.true, n.fine, n.fine )
    # Convert to matrices for plotting

image( 1:n.coarse, 1:n.coarse, m.f.w )
    # f evaluated on the approximating grid
image( xy.fine, xy.fine, m.f.v )
contour( xy.fine, xy.fine, m.f.v, add=TRUE )
    # The linear interpolation based on the coarse grid
image( xy.fine, xy.fine, m.f.true )
contour( xy.fine, xy.fine, m.f.true, add=TRUE )
    # The true function evaluated on the fine grid
}