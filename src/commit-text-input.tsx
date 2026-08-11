import { useEffect, useRef, useState, type InputHTMLAttributes, type KeyboardEvent, type TextareaHTMLAttributes } from "react";

interface CommitTextInputProps extends Omit<InputHTMLAttributes<HTMLInputElement>, "value" | "onChange"> {
  value: string;
  onCommit: (value: string) => void;
}

/** Keep keystroke state local so expensive parent visualizations do not rerender until commit. */
export function CommitTextInput({ value, onCommit, onBlur, onKeyDown, ...props }: CommitTextInputProps) {
  const [draft, setDraft] = useState(value);
  const focused = useRef(false);

  useEffect(() => {
    if (!focused.current) setDraft(value);
  }, [value]);

  return <input
    {...props}
    value={draft}
    onFocus={() => { focused.current = true; }}
    onChange={(event) => setDraft(event.target.value)}
    onBlur={(event) => {
      focused.current = false;
      onCommit(draft);
      onBlur?.(event);
    }}
    onKeyDown={(event: KeyboardEvent<HTMLInputElement>) => {
      if (event.key === "Enter") event.currentTarget.blur();
      else if (event.key === "Escape") {
        setDraft(value);
        event.currentTarget.blur();
      }
      onKeyDown?.(event);
    }}
  />;
}

interface CommitTextareaProps extends Omit<TextareaHTMLAttributes<HTMLTextAreaElement>, "value" | "onChange"> {
  value: string;
  onCommit: (value: string) => void;
}

/** Multiline counterpart: blur or Ctrl/Cmd+Enter commits; ordinary Enter inserts a line. */
export function CommitTextarea({ value, onCommit, onBlur, onKeyDown, ...props }: CommitTextareaProps) {
  const [draft, setDraft] = useState(value);
  const focused = useRef(false);

  useEffect(() => {
    if (!focused.current) setDraft(value);
  }, [value]);

  return <textarea
    {...props}
    value={draft}
    onFocus={() => { focused.current = true; }}
    onChange={(event) => setDraft(event.target.value)}
    onBlur={(event) => {
      focused.current = false;
      onCommit(draft);
      onBlur?.(event);
    }}
    onKeyDown={(event) => {
      if (event.key === "Enter" && (event.ctrlKey || event.metaKey)) event.currentTarget.blur();
      else if (event.key === "Escape") {
        setDraft(value);
        event.currentTarget.blur();
      }
      onKeyDown?.(event);
    }}
  />;
}
