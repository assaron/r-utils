#' @export
check_staged <- function(
        repo = ".",
        check_dir = file.path(tempdir(), "check_staged"),
        args = c("--no-manual", "--as-cran"),
        error_on = c("error", "warning", "note"),
        quiet = FALSE
) {
    error_on <- match.arg(error_on)
    
    if (!requireNamespace("rcmdcheck", quietly = TRUE)) {
        stop("Package 'rcmdcheck' is required. Install it first.")
    }
    
    git <- Sys.which("git")
    if (git == "") stop("Could not find 'git' on PATH. Install git and retry.")
    
    repo <- normalizePath(repo, winslash = "/", mustWork = TRUE)
    
    root <- suppressWarnings(
        system2(git, c("-C", repo, "rev-parse", "--show-toplevel"),
                stdout = TRUE, stderr = TRUE)
    )
    if (length(root) != 1 || !nzchar(root)) stop("Not a git repository: ", repo)
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    
    desc_path <- file.path(root, "DESCRIPTION")
    if (!file.exists(desc_path)) stop("No DESCRIPTION found at repo root: ", root)
    
    # Create a clean temp package tree: HEAD + staged changes
    pkg_tree <- tempfile("staged_pkg_tree_")
    dir.create(pkg_tree, recursive = TRUE)
    
    on.exit(unlink(pkg_tree, recursive = TRUE, force = TRUE), add = TRUE)
    
    # Export HEAD into pkg_tree (tracked files only)
    cmd <- sprintf(
        "%s -C %s archive --format=tar HEAD | tar -x -C %s",
        shQuote(git), shQuote(root), shQuote(pkg_tree)
    )
    rc <- system(cmd, ignore.stdout = quiet, ignore.stderr = quiet)
    if (rc != 0) {
        stop("Failed to export HEAD via 'git archive'. ",
             "This command requires a working 'tar' in your environment.")
    }
    
    # Create patch from staged changes (index vs HEAD)
    patch_lines <- system2(
        git,
        c("-C", root, "diff", "--cached", "--binary", "HEAD"),
        stdout = TRUE,
        stderr = TRUE
    )
    
    if (length(patch_lines) > 0 && any(nzchar(patch_lines))) {
        patch_file <- tempfile("staged_patch_", fileext = ".patch")
        writeLines(patch_lines, patch_file)
        on.exit(unlink(patch_file, force = TRUE), add = TRUE)
        
        rc2 <- system2(
            git,
            c("-C", pkg_tree, "apply", "--whitespace=nowarn", patch_file),
            stdout = if (quiet) FALSE else "",
            stderr = if (quiet) FALSE else ""
        )
        if (!identical(rc2, 0L)) {
            stop("Failed to apply staged patch to exported tree. ",
                 "Your index may include changes that do not apply cleanly to HEAD.")
        }
    }
    
    # Prepare check output directory
    check_dir <- normalizePath(check_dir, winslash = "/", mustWork = FALSE)
    dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Run R CMD check on the staged tree
    res <- rcmdcheck::rcmdcheck(
        path = pkg_tree,
        args = args,
        error_on = error_on,
        check_dir = check_dir,
        quiet = quiet
    )
    
    invisible(res)
}