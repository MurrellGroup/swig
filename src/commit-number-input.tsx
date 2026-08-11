import { useEffect, useRef, useState, type InputHTMLAttributes, type KeyboardEvent } from "react";

interface CommitNumberInputProps extends Omit<InputHTMLAttributes<HTMLInputElement>, "type" | "value" | "onChange"> {
  value: number;
  onCommit: (value: number) => void;
  blankWhenZero?: boolean;
}

/**
 * Numeric input with an unconstrained editing draft. This avoids converting
 * transient values such as "", "0." or "1e-" on every keystroke. The model
 * is updated only on Enter or blur; Escape restores the last applied value.
 */
export function CommitNumberInput({ value, onCommit, blankWhenZero = false, min, max, onBlur, onKeyDown, ...props }: CommitNumberInputProps) {
  const displayValue = (next: number) => blankWhenZero && next === 0 ? "" : String(next);
  const [draft, setDraft] = useState(displayValue(value));
  const focused = useRef(false);

  useEffect(() => {
    if (!focused.current) setDraft(displayValue(value));
  }, [value, blankWhenZero]);

  const commit = () => {
    const parsed = Number(draft);
    if (!draft.trim() && blankWhenZero) {
      onCommit(0);
      setDraft("");
      return;
    }
    if (!draft.trim() || !Number.isFinite(parsed)) {
      setDraft(displayValue(value));
      return;
    }
    const minimum = min === undefined ? Number.NEGATIVE_INFINITY : Number(min);
    const maximum = max === undefined ? Number.POSITIVE_INFINITY : Number(max);
    const next = Math.min(maximum, Math.max(minimum, parsed));
    onCommit(next);
    setDraft(displayValue(next));
  };

  return <input
    {...props}
    type="number"
    min={min}
    max={max}
    value={draft}
    onFocus={() => { focused.current = true; }}
    onChange={(event) => setDraft(event.target.value)}
    onBlur={(event) => {
      focused.current = false;
      commit();
      onBlur?.(event);
    }}
    onKeyDown={(event: KeyboardEvent<HTMLInputElement>) => {
      if (event.key === "Enter") event.currentTarget.blur();
      else if (event.key === "Escape") {
        setDraft(displayValue(value));
        event.currentTarget.blur();
      }
      onKeyDown?.(event);
    }}
  />;
}
