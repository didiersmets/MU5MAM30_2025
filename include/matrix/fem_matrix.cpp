#include "fem_matrix.h"
#include "sys_utils.h"

void FEMatrix::mvp(const double *__restrict x, double *__restrict y) const
{
    ASSERT(fem_type == P1_sym);
    ASSERT(m != nullptr);
    
    size_t n = m->vertex_count();
    
    // Diagonal contribution
    for (size_t i = 0; i < n; i++) {
        y[i] = diag[i] * x[i];
    }
    
    // Off-diagonal contributions
    for (size_t t = 0; t < m->triangle_count(); t++) {
        uint32_t v0 = m->indices[3*t + 0];
        uint32_t v1 = m->indices[3*t + 1];
        uint32_t v2 = m->indices[3*t + 2];
        
        size_t base = 3 * t;
        double s01 = off_diag[base + 0];
        double s02 = off_diag[base + 1];
        double s12 = off_diag[base + 2];
        
        y[v0] += s01 * x[v1] + s02 * x[v2];
        y[v1] += s01 * x[v0] + s12 * x[v2];
        y[v2] += s02 * x[v0] + s12 * x[v1];
    }
}

double FEMatrix::sum() const
{
    ASSERT(fem_type == P1_sym);
    
    double total = 0.0;
    
    for (size_t i = 0; i < diag.size; i++) {
        total += diag[i];
    }
    
    for (size_t i = 0; i < off_diag.size; i++) {
        total += 2.0 * off_diag[i];
    }
    
    return total;
}