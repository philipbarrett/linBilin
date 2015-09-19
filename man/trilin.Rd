\name{trilin}
\alias{trilin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast 3D linear interpolation
}
\description{
Takes values of a function evaluated on a grid and computes an approximate value on the interior of the grid using linear interpolation.
}
\usage{
trilin(pt, x, y, z, m_f, n_x, n_y, n_z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pt}{An \eqn{(x,y,z)} vector at which to approximately evaluate a function}
  \item{x}{A vector of x-values}
  \item{y}{A vector of y-values}
  \item{z}{A vector of z-values}
  \item{m_f}{The matrix of function values on the grid \eqn{x \times y \times z}.  Must be ordered vertically by \code{x} and then \code{y}, and horizontally by \code{z}.  For example, if \eqn{x=(x_1,x_2)}, \eqn{y=(y_1,y_2)}, and \eqn{z=(z_1,z_2)}, then \code{m_f} is:
  \deqn{R=\left[  \begin{array}{cc}
                    f(x_1, y_1, z_1)  & f(x_1, y_1, z_2) \\
                    f(x_1, y_2, z_1)  & f(x_1, y_2, z_2) \\
                    f(x_2, y_1, z_1)  & f(x_2, y_1, z_2) \\
                    f(x_2, y_2, z_1)  & f(x_2, y_2, z_2)
                  \end{array}\right]}}
  \item{n_x}{An integer.  The length of \code{x}}
  \item{n_y}{An integer.  The length of \code{y}}
  \item{n_z}{An integer.  The length of \code{z}}
}

\value{
Returns the trilinear interpolation of \eqn{(x,y,z,f)} evaluated at \code{pt}}.
}
\examples{
f <- function(pt) pt[1] + pt[2] * pt[3] - pt[3] * pt[1]
    # An arbitrary function
nn <- 7
    # Number of points at which the interpolating grid is defined
v.x <- v.y <- v.z <- 1:nn
    # The vector of x, y, and z values (need not all be the same)
m.xy <- cbind( rep( v.x, each=nn), rep( v.y, nn ) )
    # The matrix of (x,y) values at which to evaluate f
v.f <- t( apply( m.xy, 1, function( xy ) sapply( v.z, function(z) f( c( xy, z ) ) ) ) )
    # The matrix of the z-values

test.pt <- 1:3
print(f(test.pt))
print(trilin( test.pt, v.x, v.y, v.z, v.f, 7, 7, 7 ) )
    # Is exact at grid points

test.pt.2 <- 1:3 + .2
print(f(test.pt.2))
print(trilin( test.pt.2, v.x, v.y, v.z, v.f, 7, 7, 7 ) )
    # Approximation is still exact with quadratic functions

}
