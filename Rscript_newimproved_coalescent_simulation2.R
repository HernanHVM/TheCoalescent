# R Script mejorado para Simulaciones Coalescentes por HVM usando Claude Sonnet 3.7 20250328
# Este Rcript simula datos genéticos bajo el modelo de coalescencia
# Creative Commons Licence 4.0

# Instala los paquetes requeridos si no están disponibles
if (!require("ape")) install.packages("ape")
if (!require("phyclust")) install.packages("phyclust")
if (!require("coala")) install.packages("coala")
library(ape)
library(phyclust)
library(coala)
activate_ms()
#########################################################
# 1. Simulación de coalescencia simple usando ape
#########################################################

simple_coalescent <- function(n, seed = NULL) {
  # n = sample size (number of individuals)
  # Returns a coalescent tree using ape's built-in function
  
  if (!is.null(seed)) set.seed(seed)
  
  # Use ape's rcoal function to simulate a coalescent tree
  tree <- rcoal(n)
  
  # Add meaningful tip labels
  tree$tip.label <- paste("Ind", 1:n, sep="")
  
  return(tree)
}

# Graficar el árbol de coalescencia usando ape y su función de plot clásica
plot_coalescent <- function(tree) {
  plot(tree, main = "Simulated Coalescent Tree", cex = 0.8)
  return(tree)
}

#########################################################
# 2. Simulación suando el paquete coala
#########################################################

# Simulación con un tamaño poblacional constante
simulate_constant_pop <- function(sample_size = 10, mutation_rate = 5, sequence_length = 1000) {
  model <- coal_model(sample_size = sample_size, loci_number = 1, loci_length = sequence_length) +
           feat_mutation(mutation_rate) +
           sumstat_trees() +
           sumstat_nucleotide_div()
  
  result <- simulate(model)
  return(result)
}

# Simulación con crecimiento poblacional
simulate_pop_growth <- function(sample_size = 10, mutation_rate = 5, sequence_length = 1000, 
                               growth_rate = 0.1, current_size = 1) {
  model <- coal_model(sample_size = sample_size, loci_number = 1, loci_length = sequence_length) +
           feat_mutation(mutation_rate) +
           feat_growth(growth_rate, current_size) +
           sumstat_trees() +
           sumstat_nucleotide_div()
  
  result <- simulate(model)
  return(result)
}

# Simulación con estructura poblacional 
simulate_pop_structure <- function(sample_sizes = c(5, 5), mutation_rate = 5, sequence_length = 1000) {
  # Create a migration matrix with proper numeric values
  # For two populations, a 2x2 matrix where diagonal is 0
  migration_rate <- 0.5
  migration_matrix <- matrix(migration_rate, nrow = 2, ncol = 2)
  diag(migration_matrix) <- 0  # No migration within same population
  
  model <- coal_model(sample_size = sample_sizes, loci_number = 1, loci_length = sequence_length) +
           feat_mutation(mutation_rate) +
           feat_migration(0.5, symmetric = TRUE) +  #migration_matrix
           sumstat_trees()
  
  # Try to add SFS and nucleotide diversity separately
  tryCatch({
    model <- model + sumstat_nucleotide_div()
  }, error = function(e) {
    cat("Could not add nucleotide diversity:", e$message, "\n")
  })
  
  tryCatch({
    model <- model + sumstat_sfs()
  }, error = function(e) {
    cat("Could not add SFS calculation:", e$message, "\n")
  })
  
  result <- simulate(model)
  return(result)
}

# Alternative approach for structured populations
simulate_structured_alt <- function(sample_sizes = c(5, 5), 
                                   migration_rate = 0.5, 
                                   mutation_rate = 5,
                                   sequence_length = 1000) {
  
  # Number of populations
  n_pops <- length(sample_sizes)
  
  # Try a simpler approach with island model
  tryCatch({
    model <- coal_model(sample_size = sample_sizes, 
                       loci_number = 1, 
                       loci_length = sequence_length) +
             feat_mutation(mutation_rate) +
             feat_structure(value = n_pops, migration_rate) +
             sumstat_trees()
    
    # Try to add other statistics
    tryCatch({
      model <- model + sumstat_nucleotide_div()
    }, error = function(e) {
      cat("Could not add nucleotide diversity:", e$message, "\n")
    })
    
    tryCatch({
      model <- model + sumstat_sfs()
    }, error = function(e) {
      cat("Could not add SFS calculation:", e$message, "\n")
    })
    
    result <- simulate(model)
    return(result)
  }, error = function(e) {
    cat("Island model failed:", e$message, "\n")
    cat("Falling back to basic simulation approach...\n")
    
    # Fallback to basic approach - simulate populations separately
    return(simulate_basic_structured(sample_sizes, mutation_rate, sequence_length))
  })
}

# Basic structured population simulation as a fallback
simulate_basic_structured <- function(sample_sizes = c(5, 5), 
                                     mutation_rate = 5,
                                     sequence_length = 1000) {
  # This is a simplified approach that simulates populations separately
  
  results <- list()
  
  # Simulate each population separately
  for (i in 1:length(sample_sizes)) {
    pop_model <- coal_model(sample_size = sample_sizes[i], 
                           loci_number = 1, 
                           loci_length = sequence_length) +
                 feat_mutation(mutation_rate) +
                 sumstat_trees() +
                 sumstat_nucleotide_div()
    
    results[[i]] <- simulate(pop_model)
  }
  
  # Create a simple combined result
  combined <- list(
    pop_results = results,
    pi = mean(sapply(results, function(x) x$pi)),
    trees = lapply(results, function(x) x$trees[[1]])
  )
  
  return(combined)
}

# Función de ayuda para graficar árboles de coala 
# Esto arregla el problema con árboles de coala que no graficaban bien 
plot_coala_tree <- function(coala_result, main_title = "Árbol de Coalescencia") {
  # Extraer el árbol con el resultado de coala 
  if (is.null(coala_result$trees) || length(coala_result$trees) == 0) {
    cat("No se encontraron árboles en el resultado de la simulación\n")
    return(NULL)
  }
  
  # Convertirlo a phylo object de ape si es necesario
  coala_tree <- coala_result$trees[[1]]
  
  # Si ya es un phylo object, graficarlo
  if (inherits(coala_tree, "phylo")) {
    plot(coala_tree, main = main_title)
    return(coala_tree)
  }
  
  # Si ya tenemos una tira de texto newick, convertirla 
  if (is.character(coala_tree)) {
    tree <- read.tree(text = coala_tree)
    if (!is.null(tree)) {
      plot(tree, main = main_title)
      return(tree)
    }
  }
  
  # Si no podemos graficar directamente, generar un arbolito al azar con igual número de puntas 
  cat("Could not directly plot coala tree, generating substitute...\n")
  n_tips <- sample_size <- 10  # default
  if (!is.null(coala_result$model$sample_size)) {
    if (is.list(coala_result$model$sample_size)) {
      n_tips <- sum(unlist(coala_result$model$sample_size))
    } else {
      n_tips <- coala_result$model$sample_size
    }
  }
  
  # Generar un árbol de coalescencia al azar con el mismo número de puntas 
  tree <- rcoal(n_tips)
  plot(tree, main = paste(main_title, "(Substitute)"))
  
  return(tree)
}

#########################################################
# 3. Correr Simulaciones de Ejemplo
#########################################################

# Configurar un número semilla para reproducibilidad
set.seed(42)
  
# Restaurar panel de graficación
par(mfrow = c(1, 1))

# Ejemplo 1: Simular de coalescencia simple y visualización
cat("\nCorriendo simulación de coalescencia simple...\n")
tree <- simple_coalescent(10)
plot_coalescent(tree)

# Ejemplo 2: Simular con tamaño poblacional constante
cat("\nCorriendo simulación con tamaño poblacional constante...\n")
const_sim <- simulate_constant_pop()
cat("Diversidad nucleotídica:", const_sim$pi, "\n")
# Plot the tree from this simulation using the fixed function
plot_coala_tree(const_sim, "Árbol de Coalescencia (Población Constante)")

# Ejemplo 3: Simular con crecimiento poblacional
cat("\nCorriendo simulación con crecimiento constante...\n")
growth_sim <- simulate_pop_growth()
cat("Diversidad nucleotídica (población en crecimiento):", growth_sim$pi, "\n")
# Grafica el árbol de ésta simulación
plot_coala_tree(growth_sim, "Árbol de Coalescencia (Crecimiento Poblacional)")

# Ejemplo 4: Simular con estructura poblacional (dos poblaciones con migración)
# Tratar la versión corregida primero, usar alternativas si falla 
cat("\nRunning simulation with population structure...\n")
struct_sim <- NULL
tryCatch({
  struct_sim <- simulate_pop_structure()
}, error = function(e) {
  cat("Original structured simulation failed:", e$message, "\n")
  cat("Trying alternative implementation...\n")
  struct_sim <<- simulate_structured_alt()
})

if (!is.null(struct_sim)) {
  # Check if pi exists in the result
  if (!is.null(struct_sim$pi)) {
    cat("Nucleotide diversity (structured population):", struct_sim$pi, "\n")
  }
  
  # Check if sfs exists in the result
  if (!is.null(struct_sim$sfs)) {
    cat("SFS (Site Frequency Spectrum aka derived alleles):", struct_sim$sfs, "\n")
  }
  
  # Plot the tree from this simulation
  plot_coala_tree(struct_sim, "Coalescent Tree (Structured Population)")
} else {
  cat("Failed to run structured population simulation\n")
}

#########################################################
# 4. Compara diferentes escenarios demográficos
#########################################################

compare_scenarios <- function() {
  # Correr múltiples simulaciones y recolectar estadísticas 
  n_sims <- 10  # Reducir de 20 para ahorrar tiempo
  
  # Almacenaje de resultado
  results <- data.frame(
    Scenario = character(n_sims * 3),
    NucleotideDiversity = numeric(n_sims * 3),
    Time2MRCA = numeric(n_sims * 3)
  )
  
  # Correr simulaciones 
  for (i in 1:n_sims) {
    # Población constante 
    const_sim <- simulate_constant_pop()
    results$Scenario[(i-1)*3 + 1] <- "Constante"
    results$NucleotideDiversity[(i-1)*3 + 1] <- const_sim$pi
    
    # Usar otro acercamiento para obtener el TMRCA de forma más robusta
    tree <- NULL
    tryCatch({
      if (!is.null(const_sim$trees) && length(const_sim$trees) > 0) {
        if (inherits(const_sim$trees[[1]], "phylo")) {
          tree <- const_sim$trees[[1]]
        } else if (is.character(const_sim$trees[[1]])) {
          tree <- read.tree(text = const_sim$trees[[1]])
        }
      }
    }, error = function(e) {
      # Si no podemos extraer un árbol apropiado, usemos NA para el TMRCA 
      cat("Error processing tree:", e$message, "\n")
    })
    
    if (!is.null(tree)) {
      results$Time2MRCA[(i-1)*3 + 1] <- max(node.depth.edgelength(tree))
    } else {
      results$Time2MRCA[(i-1)*3 + 1] <- NA
    }
    
    # Población en crecimiento 
    growth_sim <- simulate_pop_growth()
    results$Scenario[(i-1)*3 + 2] <- "Creci"
    results$NucleotideDiversity[(i-1)*3 + 2] <- growth_sim$pi
    
    # Acercamiento similar al escenario de crecimiento 
    tree <- NULL
    tryCatch({
      if (!is.null(growth_sim$trees) && length(growth_sim$trees) > 0) {
        if (inherits(growth_sim$trees[[1]], "phylo")) {
          tree <- growth_sim$trees[[1]]
        } else if (is.character(growth_sim$trees[[1]])) {
          tree <- read.tree(text = growth_sim$trees[[1]])
        }
      }
    }, error = function(e) {
      cat("Error processing tree:", e$message, "\n")
    })
    
    if (!is.null(tree)) {
      results$Time2MRCA[(i-1)*3 + 2] <- max(node.depth.edgelength(tree))
    } else {
      results$Time2MRCA[(i-1)*3 + 2] <- NA
    }
    
    # Población estructurada, utilizar implementación alternativa que es más robusta 
    struct_sim <- NULL
    tryCatch({
      #struct_sim <- simulate_structured_alt() # alt implementation
      struct_sim <- simulate_pop_structure()
    }, error = function(e) {
      cat("Original structured simulation failed:", e$message, "\n")
    })
    
    results$Scenario[(i-1)*3 + 3] <- "Estructurada"
    
    if (!is.null(struct_sim) && !is.null(struct_sim$pi)) {
      results$NucleotideDiversity[(i-1)*3 + 3] <- struct_sim$pi
    } else {
      results$NucleotideDiversity[(i-1)*3 + 3] <- NA
    }
    
    # Para el escenario estructurado 
    tree <- NULL
    if (!is.null(struct_sim)) {
      tryCatch({
        if (!is.null(struct_sim$trees) && length(struct_sim$trees) > 0) {
          if (inherits(struct_sim$trees[[1]], "phylo")) {
            tree <- struct_sim$trees[[1]]
          } else if (is.character(struct_sim$trees[[1]])) {
            tree <- read.tree(text = struct_sim$trees[[1]])
          }
        }
      }, error = function(e) {
        cat("Error processing tree:", e$message, "\n")
      })
    }
    
    if (!is.null(tree)) {
      results$Time2MRCA[(i-1)*3 + 3] <- max(node.depth.edgelength(tree))
    } else {
      results$Time2MRCA[(i-1)*3 + 3] <- NA
    }
  }
  
  # Remover NAs para graficar 
  results_clean <- results[!is.na(results$Time2MRCA) & !is.na(results$NucleotideDiversity),]
  
  # Crear boxplots para comparar escenarios si tenemos suficientes datos 
  if (nrow(results_clean) > 0) {
    # Revisar si tenemos los tres escenarios 
    scenarios <- unique(results_clean$Scenario)
    if (length(scenarios) >= 2) {
      par(mfrow=c(1,2))
      boxplot(NucleotideDiversity ~ Scenario, data=results_clean, 
              main="Diversidad Nucleotidica X Escenario",
              col=c("blue", "green", "orange"))
      
      boxplot(Time2MRCA ~ Scenario, data=results_clean, 
              main="Tiempo al MRCA",
              col=c("blue", "green", "orange"))
    } else {
      cat("Not enough different scenarios with valid data to create comparative boxplots\n")
    }
  } else {
    cat("Not enough valid data to create boxplots\n")
  }
  
  return(results)
}

# Correr la comparación 
cat("\nComparando diferentes escenarios demograficos...\n")
comparison_results <- compare_scenarios()

# Imprimir estadísticas de resumen 
cat("\nResumen de los resultados comparados:\n")
print(tapply(comparison_results$NucleotideDiversity, comparison_results$Scenario, mean, na.rm=TRUE))
print(tapply(comparison_results$Time2MRCA, comparison_results$Scenario, mean, na.rm=TRUE))

# Regresar el panel de graficación al original 1x1 
par(mfrow=c(1,1))

#########################################################
# 5. Funciones adicionales para explorar propiedades de la coalescencia 
#########################################################

# Simular y visualizar el efecto del tamaño de muestra en los árboles de coalescencia 
explore_sample_size <- function(sample_sizes = c(5, 10, 20, 50)) {
  # Partir panel de graficacion 2x2 
  par(mfrow = c(2, 2))
  
  # Simular y graficar árboles de diferentes tamaños de muestra 
  for (n in sample_sizes) {
    tree <- simple_coalescent(n)
    plot(tree, main = paste("n =", n), cex = 0.6)
    
    # Agregar una barra de escala temporal
    h <- max(node.depth.edgelength(tree))
    cat("Altura del arbol para n =", n, ":", h, "\n")
  }
  
  # Regresar el panel de graficación al original 1x1 
  par(mfrow = c(1, 1))
}

# Investigar la distribución del tiempo al MRCA (TMRCA). 
investigate_tmrca <- function(n_sims = 100, sample_size = 10) {
  # Almacenar los valores de TMRCA
  tmrca_values <- numeric(n_sims)
  
  # Correr simulaciones 
  for (i in 1:n_sims) {
    tree <- simple_coalescent(sample_size)
    tmrca_values[i] <- max(node.depth.edgelength(tree))
  }
  
  # Graficar histograma de valores de TMRCA
  hist(tmrca_values, breaks = 20, main = "Distribucion de TMRCA",
       xlab = "Tiempo al MRCA", col = "lightblue")
  
  # Agregar la expectativa teórica 
  # Para una coalescencia estándar, E[TMRCA] = 2(1-1/n)
  expected_tmrca <- 2 * (1 - 1/sample_size)
  abline(v = expected_tmrca, col = "red", lwd = 2)
  text(expected_tmrca*1.1, max(hist(tmrca_values, breaks = 20, plot = FALSE)$counts)*0.9, 
       paste("Esperado =", round(expected_tmrca, 2)), col = "red")
  
  # Mostrar estadísticas de resumen 
  cat("Resumen de valores TMRCA:\n")
  print(summary(tmrca_values))
  cat("TMRCA esperados:", expected_tmrca, "\n")
  
  return(tmrca_values)
}

# Mostrar estas funciones adicionales 
cat("\nExplorando el efecto del tamaño de muestra en árboles de coalescencia...\n")
explore_sample_size()

# Regresar el panel de graficación al original 1x1 
par(mfrow = c(1, 1))

cat("\nInvestigando la distribucion de TMRCA...\n")
tmrca_values <- investigate_tmrca()

cat("\n¡Simulacion completa!\n")
