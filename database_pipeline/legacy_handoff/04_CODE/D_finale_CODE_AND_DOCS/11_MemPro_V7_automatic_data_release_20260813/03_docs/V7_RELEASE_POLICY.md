# MemPro V7 automatic data release policy

Human manual review was omitted by explicit user approval. Records enter PUBLIC only through deterministic identity, foreign-key, evidence and business-rule checks. Records that cannot be uniquely assigned are not guessed: they are retained in FROZEN_REVIEW or EXCLUDED_AUDIT. A/B/C precedence is A > B > C while secondary mechanisms are retained. Positive and negative evidence never overwrite one another. Missing values are never converted to zero.
