# Shared logger setup for the TeoNAM 118K pipeline scripts.
# Source once, right after ROOT/setwd:  source(file.path(ROOT, "scripts/logging.R"))
#
# Emits timestamped, sprintf-style logs:  [HH:MM:SS] INFO: message
# Use log_info() / log_warn() / log_error() with "%s"/"%d"/"%.1f" formats (NOT
# paste/sprintf inside -- the formatter is sprintf). For any family/coverage/lambda
# loop, log a running ETA after each iteration, e.g.:
#   log_info(">>> %d/%d done | elapsed %.1f min | avg %.1f min | ETA ~%.1f min remaining",
#            i, N, el, el / i, (el / i) * (N - i))
# where el <- as.numeric(difftime(Sys.time(), t0, units = "mins")).
suppressMessages(library(logger))
log_layout(layout_glue_generator(format = '[{format(time, "%H:%M:%S")}] {level}: {msg}'))
log_formatter(formatter_sprintf)
log_threshold(INFO)
