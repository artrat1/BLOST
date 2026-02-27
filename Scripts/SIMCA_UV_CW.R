# ============================================================
# POPRAWIONA ADAPTACYJNA CV - WYRÓWNANIE DŁUGOŚCI WEKTORÓW
# ============================================================

library(mdatools)
library(caret)

# ============================================================
# 1. WCZYTANIE DANYCH
# ============================================================
sink("wynik.txt")
try(setwd ("c:/Users/ar/dropbox/Kamila/paper_2025/ResultsChem/Revision/Kod/SIMCA_UV_CV/")) # Dom AR
try(setwd ("c:/Users/ar123/dropbox/Kamila/paper_2025/ResultsChem/Revision/Kod/SIMCA_UV_CV/"))  
df1=read.csv2("dataUV.csv", head=FALSE, dec = ",")
data <- read.csv2("dataUV.csv", header = TRUE, sep = ";", dec = ",", 
                  stringsAsFactors = FALSE, check.names = FALSE)
data_calib <- data[data$C2 == "kalib", ]

X <- as.matrix(data_calib[, 4:ncol(data_calib)])
y <- data_calib$plant

cat("=== DANE ===\n")
print(table(y))
cat("\n")

# ============================================================
# 2. PARAMETRY CV
# ============================================================

n_repeats <- 5
classes <- c("L", "S", "B", "O", "T")

# Konfiguracja dla każdej klasy
class_config <- list(
  "L" = list(ncomp = 3, cv_method = "cv", cv_folds = 7),
  "S" = list(ncomp = 5, cv_method = "cv", cv_folds = 10),
  "B" = list(ncomp = 4, cv_method = "cv", cv_folds = 7),
  "O" = list(ncomp = 4, cv_method = "cv", cv_folds = 7),
  "T" = list(ncomp = 2, cv_method = "loo", cv_folds = NA)
)

# ============================================================
# 3. GŁÓWNA PĘTLA CV
# ============================================================

all_results <- list()

for(rep in 1:n_repeats) {
  
  set.seed(123 + rep)
  cat(sprintf("\n--- Powtórzenie %d/%d ---\n", rep, n_repeats))
  
  # Inicjalizacja pustych wektorów
  rep_pred <- character(0)
  rep_true <- character(0)
  
  # Dla klas z CV - tworzymy foldy
  n_folds <- 10  # Maksymalna liczba foldów
  fold_predictions <- list()
  fold_truth <- list()
  
  for(fold in 1:n_folds) {
    
    # Indeksy testowe dla klas innych niż T (CV)
    non_T_idx <- which(y != "T")
    non_T_y <- y[non_T_idx]
    
    # Tworzymy foldy tylko raz na powtórzenie
    if(fold == 1) {
      cv_folds_non_T <- createFolds(non_T_y, k = 10, list = TRUE, returnTrain = FALSE)
    }
    
    # Dla klasy T - LOO (rotacyjnie)
    T_idx <- which(y == "T")
    n_T <- length(T_idx)
    test_T <- T_idx[(fold - 1) %% n_T + 1]
    
    # Indeksy testowe dla pozostałych klas
    test_non_T <- non_T_idx[cv_folds_non_T[[fold]]]
    
    # Połączenie indeksów testowych
    test_idx <- c(test_non_T, test_T)
    train_idx <- setdiff(1:nrow(X), test_idx)
    
    # Sprawdzenie czy w treningowym są wszystkie klasy
    y_train <- y[train_idx]
    missing_classes <- setdiff(classes, unique(y_train))
    
    if(length(missing_classes) > 0) {
      cat(sprintf("    Fold %d: brak klas %s - pomijam\n", 
                  fold, paste(missing_classes, collapse = ", ")))
      next
    }
    
    X_train <- X[train_idx, ]
    X_test <- X[test_idx, ]
    y_test <- y[test_idx]
    
    # Trenowanie modeli
    models <- list()
    
    for(cls in classes) {
      X_cls <- X_train[y_train == cls, , drop = FALSE]
      
      # Bezpieczna liczba komponentów
      ncomp <- min(class_config[[cls]]$ncomp, nrow(X_cls) - 1)
      if(ncomp < 1) ncomp <- 1
      
      tryCatch({
        m <- simca(X_cls, cls, ncomp = ncomp)
        models[[cls]] <- selectCompNum(m, ncomp = ncomp)
      }, error = function(e) {
        print ("Error')")
      })
    }
    
    # Sprawdzenie czy udało się wytrenować wszystkie modele
    if(length(models) == length(classes)) {
      tryCatch({
        mm <- simcam(models)
        pred <- predict(mm, X_test)
        
        # Dodanie wyników tylko jeśli długości są zgodne
        if(length(pred$c.pred) == length(y_test)) {
          rep_pred <- c(rep_pred, as.character(pred$c.pred))
          rep_true <- c(rep_true, as.character(y_test))
          cat(sprintf("    Fold %d: OK (%d próbek)\n", fold, length(y_test)))
        }
      }, error = function(e) {
        cat(sprintf("    Fold %d: błąd predykcji - pomijam\n", fold))
      })
    } else {
      cat(sprintf("    Fold %d: tylko %d/%d modeli - pomijam\n", 
                  fold, length(models), length(classes)))
    }
  }
  
  # Sprawdzenie czy mamy wyniki
  if(length(rep_pred) == 0 || length(rep_true) == 0) {
    cat("  BRAK WYNIKÓW dla tego powtórzenia!\n")
    next
  }
  
  if(length(rep_pred) != length(rep_true)) {
    cat(sprintf("  Uwaga: różne długości! pred=%d, true=%d\n", 
                length(rep_pred), length(rep_true)))
    # Wyrównanie długości (obcięcie do krótszego)
    min_len <- min(length(rep_pred), length(rep_true))
    rep_pred <- rep_pred[1:min_len]
    rep_true <- rep_true[1:min_len]
  }
  
  cat(sprintf("  Łącznie: %d próbek testowych\n", length(rep_true)))
  
  # Macierz pomyłek
  tryCatch({
    cm <- confusionMatrix(
      factor(rep_pred, levels = classes),
      factor(rep_true, levels = classes)
    )
    
    all_results[[rep]] <- list(
      accuracy = cm$overall["Accuracy"],
      balanced_accuracy = mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE),
      class_metrics = data.frame(
        Class = rownames(cm$byClass),
        BalancedAccuracy = cm$byClass[, "Balanced Accuracy"],
        Precision = cm$byClass[, "Pos Pred Value"],
        Recall = cm$byClass[, "Sensitivity"],
        F1 = cm$byClass[, "F1"]
      ),
      confusion_matrix = cm$table
    )
    
    cat(sprintf("  Accuracy: %.3f\n", cm$overall["Accuracy"]))
    
  }, error = function(e) {
    cat("  Błąd tworzenia macierzy pomyłek:", e$message, "\n")
  })
}

# ============================================================
# 4. PODSUMOWANIE (tylko jeśli są wyniki)
# ============================================================

if(length(all_results) > 0) {
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("=== WYNIKI ADAPTACYJNEJ CV ===\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  cat("Konfiguracja:\n")
  cat("  - Klasa L (9 próbek): 7-fold CV\n")
  cat("  - Klasa S (24 próbki): 10-fold CV\n")
  cat("  - Klasa B (12 próbek): 7-fold CV\n")
  cat("  - Klasa O (12 próbek): 7-fold CV\n")
  cat("  - Klasa T (6 próbek): Leave-One-Out CV\n\n")
  
  cat(sprintf("Liczba udanych powtórzeń: %d/%d\n\n", length(all_results), n_repeats))
  
  # Ogólna dokładność
  accuracies <- sapply(all_results, function(x) x$accuracy)
  cat(sprintf("Accuracy: %.3f ± %.3f\n", mean(accuracies), sd(accuracies)))
  cat(sprintf("95%% CI: [%.3f, %.3f]\n\n",
              mean(accuracies) - 1.96 * sd(accuracies)/sqrt(length(accuracies)),
              mean(accuracies) + 1.96 * sd(accuracies)/sqrt(length(accuracies))))
  
  # Metryki per klasa
  cat("Metryki per klasa (średnia ± SD):\n")
  cat(rep("-", 80), "\n")
  cat(sprintf("%-8s | %-20s | %-15s | %-15s | %-10s\n",
              "Klasa", "Balanced Accuracy", "Precision", "Recall", "F1"))
  cat(rep("-", 80), "\n")
  
  for(cls in classes) {
    ba <- sapply(all_results, function(x) {
      val <- x$class_metrics$BalancedAccuracy[x$class_metrics$Class == cls]
      if(length(val) == 0) NA else val
    })
    prec <- sapply(all_results, function(x) {
      val <- x$class_metrics$Precision[x$class_metrics$Class == cls]
      if(length(val) == 0) NA else val
    })
    rec <- sapply(all_results, function(x) {
      val <- x$class_metrics$Recall[x$class_metrics$Class == cls]
      if(length(val) == 0) NA else val
    })
    f1 <- sapply(all_results, function(x) {
      val <- x$class_metrics$F1[x$class_metrics$Class == cls]
      if(length(val) == 0) NA else val
    })
    
    # Usuwanie NA
    ba <- ba[!is.na(ba)]
    prec <- prec[!is.na(prec)]
    rec <- rec[!is.na(rec)]
    f1 <- f1[!is.na(f1)]
    
    if(length(ba) > 0) {
      cat(sprintf("%-8s | %.3f ± %.3f | %.3f ± %.3f | %.3f ± %.3f | %.3f ± %.3f\n",
                  cls, 
                  mean(ba), sd(ba),
                  mean(prec), sd(prec),
                  mean(rec), sd(rec),
                  mean(f1), sd(f1)))
    } else {
      cat(sprintf("%-8s | brak danych | brak danych | brak danych | brak danych\n", cls))
    }
  }
  
  # Uśredniona macierz pomyłek
  cat("\n\n=== UŚREDNIONA MACIERZ POMYŁEK ===\n")
  avg_cm <- Reduce(`+`, lapply(all_results, function(x) x$confusion_matrix)) / length(all_results)
  print(round(avg_cm, 1))
  
} else {
  cat("\n!!! BRAK WYNIKÓW - żadne powtórzenie nie zakończyło się sukcesem !!!\n")
}

cat("\n=== KONIEC ===\n")