while IFS= read -r line; do
  if [[ "$line" == *"statistics: AS"* ]]; then
    printf '%s\n' "$line"
  fi
done < $1
