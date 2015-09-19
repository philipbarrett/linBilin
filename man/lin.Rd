\name{lin}
\alias{lin}
\title{Fast 1D Linear interpolation}
\usage{
lin( x, y, pt, n_x)
}
\arguments{
  \item{x}{A vector of x-values}
  \item{y}{A vector of y-values}
  \item{pt}{An x-value}
  \item{n_x}{An integer.  The length of \code{x} (and \code{y}).}
}
\value{
Returns the linear interpolation of (x,x) evaluated at x=pt.  That is, if \code{x}\eqn{= (x_1, x_2, \ldots, x_n)}, \code{y}\eqn{= (y_1, y_2, \ldots, y_n)}, and \code{pt}\eqn{= x,  x_i < x \le x_{i+1}}, then the output is:
\deqn{y = \left( \frac{x - x_i}{x_{i+1} - x_i} \right) y_i + \left( 1 - \frac{x - x_i}{x_{i+1} - x_i} \right) y_{i+1} }
}
\description{
  Provides very fast 1D linear interpolation.  Typically about 10 times faste than R's built-in \code{approx} function.
}
\examples{
## 1D Interpolation
f.1 <- function(x) x^2
n <- 5
x.grid <- seq( 0, 10, length.out = n )
f.1.grid <- sapply( x.grid, f.1 )
lin( x.grid, f.1.grid, 3, n ) == 3
    # Check interpolation
x.grid.fine <- seq( 0, 10, length.out = n^2 )
plot( x.grid.fine, sapply( x.grid.fine, f.1 ),
      main='Exact and linear interpoliation of x^2', type='l' )
lines( x.grid.fine, sapply( x.grid.fine, lin, x=x.grid, y=f.1.grid, n_x=n ), 
        col='red' )
    # Compare the exact and approximate solution
library(microbenchmark)
microbenchmark(sapply( x.grid.fine, lin, x=x.grid, y=f.1.grid, n_x=n ))
microbenchmark(sapply( x.grid.fine, approx, x=x.grid, y=f.1.grid ))
    # lin is about 10x faster than R's built-in approx function
}