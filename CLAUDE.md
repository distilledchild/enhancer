# Agent Instructions

> This file is mirrored across CLAUDE.md, AGENTS.md, and GEMINI.md so the same instructions load in any AI environment.

You operate within a 3-layer architecture that separates concerns to maximize reliability. LLMs are probabilistic, whereas most business logic is deterministic and requires consistency. This system fixes that mismatch.

## The 3-Layer Architecture

**Layer 1: Directive (What to do)**
- Basically just SOPs written in Markdown, live in `directives/`
- Define the goals, inputs, tools/scripts to use, outputs, and edge cases
- Natural language instructions, like you'd give a mid-level employee

**Layer 2: Orchestration (Decision making)**
- This is you. Your job: intelligent routing.
- Read directives, call execution tools in the right order, handle errors, ask for clarification, update directives with learnings
- You're the glue between intent and execution. E.g you don't try scraping websites yourself—you read `directives/scrape_website.md` and come up with inputs/outputs and then run `execution/scrape_single_site.py`

**Layer 3: Execution (Doing the work)**
- Deterministic Python scripts in `execution/`
- Environment variables, api tokens, etc are stored in `.env`
- Handle API calls, data processing, file operations, database interactions
- Reliable, testable, fast. Use scripts instead of manual work. Commented well.

**Why this works:** if you do everything yourself, errors compound. 90% accuracy per step = 59% success over 5 steps. The solution is push complexity into deterministic code. That way you just focus on decision-making.

## Operating Principles

**1. Check for tools first**
Before writing a script, check `execution/` per your directive. Only create new scripts if none exist.

**2. Self-anneal when things break**
- Read error message and stack trace
- Fix the script and test it again (unless it uses paid tokens/credits/etc—in which case you check w user first)
- Update the directive with what you learned (API limits, timing, edge cases)
- Example: you hit an API rate limit → you then look into API → find a batch endpoint that would fix → rewrite script to accommodate → test → update directive.

**3. Update directives as you learn**
Directives are living documents. When you discover API constraints, better approaches, common errors, or timing expectations—update the directive. But don't create or overwrite directives without asking unless explicitly told to. Directives are your instruction set and must be preserved (and improved upon over time, not extemporaneously used and then discarded).

## Self-annealing loop

Errors are learning opportunities. When something breaks:
1. Fix it
2. Update the tool
3. Test tool, make sure it works
4. Update directive to include new flow
5. System is now stronger

## File Organization

**Deliverables vs Intermediates:**
- **Deliverables**: Google Sheets, Google Slides, or other cloud-based outputs that the user can access
- **Intermediates**: Temporary files needed during processing

**Directory structure:**
- `.tmp/` - All intermediate files (dossiers, scraped data, temp exports). Never commit, always regenerated.
- `execution/` - Python scripts (the deterministic tools)
- `directives/` - SOPs in Markdown (the instruction set)
- `.env` - Environment variables and API keys
- `credentials.json`, `token.json` - Google OAuth credentials (required files, in `.gitignore`)

**Key principle:** Local files are only for processing. Deliverables live in cloud services (Google Sheets, Slides, etc.) where the user can access them. Everything in `.tmp/` can be deleted and regenerated.

## Summary

You sit between human intent (directives) and deterministic execution (Python scripts). Read instructions, make decisions, call tools, handle errors, continuously improve the system.

Be pragmatic. Be reliable. Self-anneal.

## Git Commit Message Conventions

When writing commit messages in this workspace:

- Do not use Conventional Commit prefixes such as `feat:`, `fix:`, `docs:`, `chore:`, or `refactor:`.
- Use a concise imperative subject line without a prefix.
- Use the commit body to describe grouped changes by workflow step when helpful.
- Keep the message factual and specific to the files staged for that commit.

## R Coding Conventions

When writing or refactoring R scripts in this workspace, follow these formatting and stylistic rules strictly:

1. **Section Headers (Step Classification)**
   Each major step (e.g., `# 1. Loop load and anchor parsing` up to `# 9. Save results`) must be enclosed above and below by a separator line of exactly 50 hash tags (`#`) with no spaces:
   ```R
   ####################################################
   # 1. Loop load and anchor parsing
   ####################################################
   ```
   Do not use spaces inside the separator line (e.g., `# ####################################################` is prohibited).

2. **Library Management & Tidyverse Single System**
   - Use `tidyverse` as the single system/framework for data analysis.
   - If other libraries (e.g., `GenomicRanges`, `data.table`) are required, they must have a comment immediately next to the `library()` call specifying the exact reason.
     Example: `library(GenomicRanges) # for genomic interval overlap and ranges operation`

3. **No `cat()` outputs, Use `print()` instead**
   - Never use `cat()` for logging or outputting results. Use `print()` exclusively.

4. **Result Annotations on print() Lines**
   - Every line containing a `print()` call must have a brief comment at the end of the line indicating ONLY the actual output value obtained from execution (without the prefix labels like "loops:" or "mean:").
     Example: `print(str_glue("loops: {nrow(loops.raw)}")) # 15085`

5. **No Solo Text print() Calls**
   - Do not use `print()` solely for outputting text labels (e.g., `print("WHERE distribution:")` is prohibited).
   - Combine text labels with the corresponding data output inside a single print statement (using `str_glue()`, `str_c()`, or string manipulation).
     Example: `print(str_glue("WHERE distribution: {str_c(str_c(names(t), t, sep=': '), collapse=', ')}"))`

6. **Tidyverse Pipe (`%>%`) Formatting**
   - If a pipe process ends with a single `%>%` operation, it can remain on one line.
   - If multiple `%>%` operations are chained (multiple steps), each step must be written on a new line with indentation.
     Example:
     ```R
     loops.w %>% 
       filter(WHERE == "UP") %>% 
       select(chr = chr1, start = start1, end = end1, loop.id, category)
     ```

7. **Comment Formatting**
   - Use only `#` for comments. Do not use dash or equality signs (e.g., `# ---`, `# ===` are prohibited).
   - For subtitle or title-like comments, use `### Title` instead of dashes.

8. **Variable and Column Naming Conventions**
    - General variables must use dot (`.`) notation in the format of `[data_type].[variable_name]` (e.g., `df.promoter` and `df.enhancer` for data frames, `gr.promoter` and `gr.enhancer` for GRanges). Avoid naming patterns like `enhancer.df` or `promoter.df` where the data type is appended at the end. Always write in the format of `[data_type].[variable_name]` (e.g. `df.enhancer` instead of `enhancer.df`).
    - File path variables must use the format `path.[extension].[description]` (e.g., `path.rds.direction` or `path.csv.loops` instead of a simple `path.rds` or `path.loop`) to clearly identify the data type and purpose.
    - Data frame column names must use underscore (`_`) notation (e.g., `test_column`, `loop_id`).

9. **File Header Description**
   - R files must start with a short header describing the file's purpose (around 5-6 lines total).
   - This block must be enclosed above and below by a separator line of exactly 50 hash tags (`#`) with no trailing characters:
     ```R
     ####################################################
     # ATAC-seq Validation of P-E Loop Anchors
     # Yuan et al. 2021 (rn7 liftover) x Hi-C loop anchors
     # Q: Is enhancer anchor in open chromatin region?
     ####################################################
     ```

10. **Working Directory Configuration (`setwd`) & Printing Options**
    - Under the library loading block, specify the R output options to prevent scientific notation:
      `options(scipen = 999) # prevent scientific notation in base R`
      `options(pillar.sigfig = 10) # display up to 10 significant digits in tibbles`
    - Immediately below, specify the project root directory explicitly using `setwd()`.
    - Provide a comment on the same line indicating the OS environment (e.g., `# Mac`, `# Ubuntu`).
    - Below `setwd()`, write `getwd()` and append a comment indicating the exact verified working directory path.
      Example:
      ```R
      library(tidyverse)
      library(GenomicRanges) # for genomic interval overlap and ranges operation

      options(scipen = 999) # prevent scientific notation in base R
      options(pillar.sigfig = 10) # display up to 10 significant digits in tibbles

      setwd("/Users/pete/Desktop/playground/enhancer") # Mac
      getwd() # /Users/pete/Desktop/playground/enhancer
      ```

11. **Function Header Documentation**
    - Every custom function must have a single-line comment directly above its definition.
    - The format must be: `# func[N]: brief English description of what the function does`, where `[N]` is the sequential number of the function defined in the file from the top (starting from 1).
    - **Crucial Rule**: The description must be written briefly in English (설명은 반드시 영문으로 간단히 작성할 것).
      Example:
      ```R
      # func1: Parse loop_id string into genomic coordinates data frame
      parse_loop_id <- function(id) { ... }

      # func2: Load narrowPeak ATAC-seq file and convert to GRanges
      load_atac <- function(path, label) { ... }
      ```

12. **Avoid `paste()` and `paste0()`**
    - Unless absolutely necessary, avoid using basic R `paste()` or `paste0()` for string concatenation.
    - Use `str_c()` or `str_glue()` from the `stringr` package instead to maintain consistency within the tidyverse environment.
      Example: Use `str_c("loops: ", nrow(df.loops))` instead of `paste("loops:", nrow(df.loops))` (note that `str_c()` does not insert spaces automatically like `paste()`).
